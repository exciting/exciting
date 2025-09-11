module calculate_dielectric_function
  use asserts, only: assert
  use constants, only : zone, zi, zzero
  use device_linalg_common_interface, only: zgemm_gpu
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo_progress, write_to_gwinfo_progress_bar, &
    progress_polarizability_tetrahedron, progress_epsilon_loop_kpoints, progress_epsilon_loop_blocks
  use gw_io, only: build_file_name, write_to_file
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
  use m_memory_device, only: allocate_device_memory, deallocate_device_memory, &
    bytes_double_complex, bytes_int, get_device_pointer
  use mod_APW_LO, only: apwordmax
  use mod_atoms, only: natmtot
  use mod_bands, only: eveck, eveckp, eveckalm, eveckpalm, &
    n_occupied_bands => nomax
  use mod_core_states, only: n_core_states => ncg
  use mod_device_offload,    only: device_world
  use mod_dielectric_function, only: epsilon, epsh, epsw1, epsw2, pmatcv, pmatvv
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_expand_products, only: expand_products_generic
  use mod_head_and_wings, only: calcwings, calchead
  use mod_misc_gw, only: Gamma
  use mod_mpi_gw, only : mpi_sum_array, indexes_parallelization
  use mod_muffin_tin, only: lmmaxapw
  use mod_pointer_remapping, only: remap_fortran_pointer
  use mod_product_basis, only: mbsiz, minmmat
  use modinput, only: input
  use modgw, only:  time_df, time_dfbody, time_dfhead, time_dfwing, &
    fnm, fnm_sum, fnm_tet_global => fnm_tet, mblksiz, kqset, Gkqset, freq
  use modmpi, only: terminate, rank
  use modxs, only : symt2
  use precision, only: dp, i32, max_length => str_64

#include "offload.fpp"

  implicit none

  private

  character(len=*), parameter :: file_name_polarizability_factor = 'POLARIZABILITY_FACTOR_Q'

  type epsilon_indexes
    type(indexes_parallelization) :: k_points
    type(indexes_parallelization) :: empty_bands
    type(indexes_parallelization) :: frequencies
  end type

  public :: calcepsilon, epsilon_indexes

contains

!> Compute the RPA dielectric matrix
subroutine calcepsilon( iq, indexes, write_progress_to_gw_info, mpi_env, print_Polarizability_Factor, file_format )
  !> Index of the q-point treated here
  integer(i32), intent(in) :: iq
  !> Set of indexes (k-point, empty bands, frequencies) used to calculate epsilon
  type(epsilon_indexes), intent(in) :: indexes
  !> If `true`, this MPI rank will write the progress of current calculation to `GW_INFO.OUT`
  logical, intent(in), optional :: write_progress_to_gw_info
  !> MPI environment
  type(mpiinfo), intent(in), optional :: mpi_env
  !> If true, print the polarizability factor to an output file.
  !> `mpi_env` must also be present to print it
  logical, intent(in), optional :: print_Polarizability_Factor
  !> Format of input/output files
  character(len=*), intent(in), optional :: file_format
  
  integer(i32) :: ie1, ie2, first_unoccupied_band, last_unoccupied_band
  integer(i32) :: iom, iomstart, iomend, ibasis
  integer(i32) :: ik, jk, ki, kf
  integer(i32) :: im, iop, jop
  integer(i32) :: ndim, mdim, nmdim
  integer(i32) :: nblk, iblk, mstart, mend
  integer(i32) :: first_occupied_band
  integer(i32) :: last_occupied_band
  real(dp)    :: tstart, tend, ti, tf, ta, tb, t_accumulate, fraction
  real(dp)    :: ti_block, tf_block, t_accumulate_block
  real(dp)    :: wto, wlo
  complex(dp) :: head(3,3), f, w
  complex(dp), pointer, contiguous :: minm(:,:,:)
  type(c_ptr) :: minm_cptr
  complex(dp), allocatable :: evecfv(:,:)
  complex(dp), allocatable :: fnm_buffer(:, :, :, :)
  complex(dp), allocatable, target :: fnm_tet(:, :, :, :)
  integer(i32) :: my_device
  logical :: print_Polarizability
  logical :: tetrahedron_method
  logical :: write_progress, write_progress_block_multiplication

  call timesec(tstart)

  call assert( present(print_Polarizability_Factor) .eqv. present(file_format), &
    'file_format and print_Polarizability_Factor must be either both present or not')

  !=============================
  ! Initialization
  !=============================
  ! Get device id
  my_device = device_world%get_device()
  fnm => null()
  ki = indexes%k_points%my_first
  kf = indexes%k_points%my_last
  first_unoccupied_band = indexes%empty_bands%my_first
  last_unoccupied_band = indexes%empty_bands%my_last
  iomstart = indexes%frequencies%my_first
  iomend = indexes%frequencies%my_last
  print_Polarizability = .false.
  if( present(mpi_env) .and. present(print_Polarizability_Factor) ) &
    & print_Polarizability = print_Polarizability_Factor
  write_progress = .false.
  if( present( write_progress_to_gw_info ) ) write_progress = write_progress_to_gw_info

  ! total number of states including the core ones
  first_occupied_band = 1
  last_occupied_band = n_occupied_bands
  if (input%gw%coreflag=='all') then
    ndim = n_occupied_bands + n_core_states
  else
    ndim = n_occupied_bands
  end if
  mdim = last_unoccupied_band - first_unoccupied_band + 1
  nmdim = ndim*mdim

  ! block size
  if (mblksiz >= mdim) then
    nblk = 1
  else
    nblk = mdim / mblksiz
    if ( mod(mdim,mblksiz) /= 0 ) nblk = nblk+1
  end if
  write_progress_block_multiplication = write_progress .and. ( nblk>1 )

  ! arrays to store products of KS eigenvectors with the matching coefficients
  allocate( eveckalm(first_occupied_band:last_occupied_band, apwordmax,lmmaxapw,natmtot) )
  allocate( eveckpalm(first_unoccupied_band:last_unoccupied_band, apwordmax,lmmaxapw,natmtot) )
  allocate( eveck(nmatmax, first_occupied_band:last_occupied_band) )
  allocate( eveckp(nmatmax, first_unoccupied_band:last_unoccupied_band) )
  allocate( evecfv(nmatmax,nstfv) )
  
  OMP_OFFLOAD target enter data map(alloc: eveckalm, eveckpalm, eveck, eveckp)

  !==================================================
  ! Calculate the q-dependent BZ integration weights
  !==================================================
  select case (trim(input%gw%qdepw))
    case('sum')
      tetrahedron_method = .false.
      allocate( fnm_sum(ndim, first_unoccupied_band:last_unoccupied_band, iomstart:iomend), source=zzero )
    case('tet')
      if( write_progress ) call write_to_gwinfo_progress( progress_polarizability_tetrahedron )
      tetrahedron_method = .true.
      allocate( fnm_tet_global(ndim, first_unoccupied_band:last_unoccupied_band, iomstart:iomend, kqset%nkpt), source=zzero )
      call qdepwtet(iq, iomstart, iomend, ndim)
      allocate( fnm_tet(ndim, first_unoccupied_band:last_unoccupied_band, iomstart:iomend, ki:kf), &
        source=fnm_tet_global(:, :, :, ki:kf) )
      deallocate( fnm_tet_global )
    case default
      call terminate( "Error(calcepsilon): Unknown qdepw method!" )
  end select
  
  if( print_Polarizability ) then
    allocate( fnm_buffer(ndim, &
                      indexes%empty_bands%global_first:indexes%empty_bands%global_last, &
                      indexes%frequencies%global_first:indexes%frequencies%global_last, &
                      indexes%k_points%global_first:indexes%k_points%global_last), source=zzero )
  end if

  !=================
  ! BZ integration
  !=================
  if( write_progress ) then 
    call timesec( tf )
    t_accumulate = 0.0_dp
    call write_to_gwinfo_progress( progress_epsilon_loop_kpoints )
  end if
  do ik = ki, kf
    if( write_progress ) ti = tf
    ! Calculate the (k,q)-dependent BZ integration weights
    if( tetrahedron_method ) then 
      fnm(1:ndim, first_unoccupied_band:last_unoccupied_band, iomstart:iomend) => fnm_tet(:, :, :, ik)
    else
      call qdepwsum(iq, ik, iomstart, iomend, ndim)
      fnm(1:ndim, first_unoccupied_band:last_unoccupied_band, iomstart:iomend) => fnm_sum
    end if
    if( print_Polarizability ) then
      fnm_buffer(:, first_unoccupied_band:last_unoccupied_band, iomstart:iomend, ik) = fnm
    end if

    !-------------------
    ! Clear memory
    !-------------------
    if (Gamma) then
      call timesec(ta)
      ! read the momentum matrix elements
      call getpmatkgw(ik)
      ! and compute the head of the dielectric function
      call calchead(ik, first_unoccupied_band, last_unoccupied_band, iomstart, iomend, ndim, epsh)
      call timesec(tb)
      time_dfhead = time_dfhead + (tb - ta)
    end if
    
    ! k-q point
    jk = kqset%kqid(ik, iq)

    ! get KS eigenvectors
    call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
    eveckp = conjg( evecfv(:, first_unoccupied_band:last_unoccupied_band) )
    call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
    eveck = evecfv(:, first_occupied_band:last_occupied_band)

    ! compute products \sum_G C_{k}n * A_{lm}
    call expand_evec(ik,'t')
    call expand_evec(jk,'c')

    OMP_OFFLOAD target update to(eveck, eveckp, eveckalm, eveckpalm)

    !=================================================
    ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
    !=================================================
    if( write_progress_block_multiplication ) then
      call timesec( tf_block )
      call write_to_gwinfo_progress( progress_epsilon_loop_blocks )
      t_accumulate_block = 0.0_dp
    end if
    do iblk = 1, nblk
      if( write_progress_block_multiplication ) ti_block = tf_block
      mstart = first_unoccupied_band + (iblk-1)*mblksiz
      mend   = min(last_unoccupied_band, mstart+mblksiz-1)
      nmdim  = ndim * (mend-mstart+1)
      allocate(minmmat(mbsiz,ndim,mstart:mend))
      OMP_OFFLOAD target enter data map(alloc: minmmat)

      ! compute M^i_{nm}+M^i_{cm}
      call expand_products_generic(ik, iq, 1, n_occupied_bands, 1, ndim-n_occupied_bands, mstart, mend, 1, 0, minmmat, .true.)
      if (Gamma) then
        call timesec(ta)
        ! wings of the dielectric matrix
        OMP_OFFLOAD target update from(minmmat)
        call calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend, mbsiz, minmmat, epsw1, epsw2)
        call timesec(tb)
        time_dfwing = time_dfwing + (tb - ta)
      end if

      ! Body
      call timesec(ta)
      call allocate_device_memory(minm_cptr, mbsiz*nmdim*bytes_double_complex, my_device)
      call c_f_pointer(minm_cptr, minm, int([mbsiz,ndim,(mend-mstart+1)],kind=c_size_t))
      ! Remapping the pointer boundaries from Fortran default
      ! TODO(mrm): When supported use lower for c_f_pointer introduced in Fortran 2023
      call remap_fortran_pointer(minm, int([1, 1, mstart], kind=i32), int([mbsiz, ndim, mend], kind=i32))
      do iom = iomstart, iomend
#if defined(FLANG_OPENMP_SLICE_MAP_BUG_WORKAROUND)
        OMP_OFFLOAD target data map(to: fnm)
#else
        OMP_OFFLOAD target data map(to: fnm(:,mstart:mend,iom))
#endif
        OMP_OFFLOAD target has_device_addr(minm)
        !$omp teams distribute parallel do collapse(3) default(none) private(ie1,ie2,ibasis) &
        !$omp shared(mstart,mend,ndim,mbsiz,minm,fnm,minmmat,iom,ik)
        do ie2 = mstart, mend
          do ie1 = 1, ndim
            do ibasis = 1, mbsiz
              minm(ibasis,ie1,ie2) = fnm(ie1,ie2,iom) * &
                                      minmmat(ibasis,ie1,ie2)
            end do
          end do ! ie1
        end do ! ie2
        !$omp end teams distribute parallel do
        OMP_OFFLOAD end target
        OMP_OFFLOAD end target data
        call zgemm_gpu( 'n', 'c', mbsiz, mbsiz, nmdim, &
                          zone, minm_cptr, mbsiz, get_device_pointer(minmmat,my_device), mbsiz, &
                          zone, get_device_pointer(epsilon(1,1,iom),my_device), mbsiz, device_world)
        call device_world%synchronize()
      end do ! iom
      call timesec(tb)
      time_dfbody = time_dfbody + (tb - ta)
      nullify(minm)
      call deallocate_device_memory(minm_cptr, my_device)
      OMP_OFFLOAD target exit data map(delete: minmmat)
      deallocate(minmmat)
      if( write_progress_block_multiplication ) then
        call timesec( tf_block )
        t_accumulate_block = t_accumulate_block + (tf_block - ti_block)
        fraction = real(iblk, dp)/(nblk)
        call write_to_gwinfo_progress_bar( t_accumulate_block, fraction, identation_level=3 )
      end if
    end do ! iblk
    if( write_progress ) then
      call timesec( tf )
      t_accumulate = t_accumulate + (tf-ti)
      fraction = real(ik - ki + 1, dp)/(kf - ki + 1)
      call write_to_gwinfo_progress_bar( t_accumulate, fraction, identation_level=2 )
    end if
  end do ! ik

  if( print_Polarizability ) then
    call mpi_sum_array( fnm_buffer, mpi_env, all_reduce=.false. )
    if( mpi_env%is_root ) call write_fnm_to_file( iq, fnm_buffer, lbound(fnm_buffer), file_format )
  end if
  
  !-------------------
  ! Clear memory
  !-------------------
  if (Gamma) then
    ! deallocate the momentum matrix elements
    deallocate(pmatvv)
    if (input%gw%coreflag=='all') deallocate(pmatcv)
  end if

  OMP_OFFLOAD target update from(epsilon)

  OMP_OFFLOAD target exit data map(delete: eveck, eveckp, eveckalm, eveckpalm)
  deallocate( eveck )
  deallocate( eveckp )
  deallocate( eveckalm )
  deallocate( eveckpalm )
  if( allocated(fnm_sum) ) deallocate( fnm_sum )
  if( associated(fnm) ) nullify( fnm )

  if (Gamma) then
    ! symmetrize \eps_{00} (head)
    do iom = iomstart, iomend
      head(:,:) = epsh(:,:,iom)
      do iop = 1, 3
        do jop = 1, 3
          call symt2app(iop, jop, 1, symt2, head, epsh(iop,jop,iom))
        end do
      end do
    end do ! iom
  end if ! Gamma

  if( present( mpi_env ) ) then
    call mpi_sum_array( epsilon, mpi_env, all_reduce=.false. )
    if( Gamma ) then
      call mpi_sum_array( epsh, mpi_env, all_reduce=.false. )
      call mpi_sum_array( epsw1, mpi_env, all_reduce=.false. )
      call mpi_sum_array( epsw2, mpi_env, all_reduce=.false. )
    end if
  end if

  !=================================
  ! \epsilon = \delta_ij - \epsilon
  !=================================
  do iom = iomstart, iomend
    if (Gamma) then
      do iop = 1, 3
        epsh(iop,iop,iom) = zone - epsh(iop,iop,iom)
      end do
    end if
    epsilon(:,:,iom) = -epsilon(:,:,iom)
    do im = 1, mbsiz
      epsilon(im,im,iom) = zone + epsilon(im,im,iom)
    end do ! im
  end do ! iom

  ! Compute contributions due to polar phonons
  if (Gamma) then
    if (rank == 0) call writedielt('EPS00', iomend-iomstart+1, freq%freqs(iomstart:iomend), epsh(:,:,iomstart:iomend), 1)
    if (input%gw%eph == 'polar') then
      ! call eph_polar(iomend-iomstart+1, cmplx(0.d0,freq%freqs(iomstart:iomend),8), epsh(:,:,iomstart:iomend))
      wlo = input%gw%wlo
      wto = input%gw%wto
      do iom = iomstart, iomend
        w = zi*freq%freqs(iom)
        ! f = (wlo**2-w**2) / (wto**2-w**2)
        f = wlo**2 / wto**2
        epsh(:,:,iom) = epsh(:,:,iom)*f
      end do
      if (rank == 0) call writedielt('EPS00+LAT', iomend-iomstart+1, freq%freqs(iomstart:iomend), epsh(:,:,iomstart:iomend), 1)
    end if
  end if

  ! timing
  call timesec(tend)
  time_df = time_df+tend-tstart

end subroutine

!> (private) Write the polarizability factor into an output file
subroutine write_fnm_to_file( idx_qpoint, polarizability_factor, lbounds, file_format )
  !> q-point index
  integer(i32), intent(in)  :: idx_qpoint 
  !> lower bounds of `polarizability_factor`
  integer(i32), intent(in)  :: lbounds(4)
  !> Polarizability factor
  complex(dp), intent(in) :: polarizability_factor(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):)
  !> Format of output file
  character(len=*), intent(in) :: file_format

  character(len=max_length) :: file_name

  call build_file_name( file_name_polarizability_factor, idx_qpoint, file_name )
  call write_to_file( file_name, polarizability_factor, lbounds, file_format )

end subroutine

end module
