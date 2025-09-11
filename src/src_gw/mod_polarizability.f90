!> This module contains procedures to store and manipulate the polarizability
module mod_polarizability
    use constants, only: zzero
    use gw_io, only: build_file_name, write_to_file, read_from_file
    use modmpi, only: terminate
    use precision, only: i32, dp
#include "offload.fpp"

    implicit none

    private

    !> polarizability function P(q). In Gamma case note that this
    !> only contains the body as the wings and the head of 
    !> the polarizability is 0.
    complex(dp), protected, allocatable :: polarizability(:,:,:)

    !> The "wings" of the polarizability.
    !> In reality they are the wings but without multiplying by the Coulomb potential
    !> This is done to prevent the recomputing of the expansions coefficients for the Gamma point.
    complex(dp), protected, allocatable ::  nocoulomb_dielectric_wing1(:,:,:)
    complex(dp), protected, allocatable ::  nocoulomb_dielectric_wing2(:,:,:)

    !> files to store the polarizability
    character(len=*), parameter :: file_name_polarizability = 'POLARIZABILITY-GW_Q'
    character(len=*), parameter :: file_name_polarizability_irreducible = 'POLARIZABILITY-GW_IQ'
    character(len=*), parameter :: file_name_nocoulomb_dielectric_wing1 = 'NOCOULOMB-EPSILON-WING1'
    character(len=*), parameter :: file_name_nocoulomb_dielectric_wing2 = 'NOCOULOMB-EPSILON-WING2'

    integer(i32), parameter   :: max_string_length = 40

    public :: write_polarizability_to_file, init_polarizability, delete_polarizability, &
              read_polarizability_from_file, compute_polarizability_at_q, from_polarizability_to_epsilon

    
contains

    subroutine init_polarizability(mbsiz, iomstart, iomend, is_Gamma_point)
        use constants, only: zzero
        implicit none
        !> Mixed basis size
        integer(i32), intent(in) :: mbsiz
        !> Frequency indexes: start and end
        integer(i32), intent(in) :: iomstart, iomend
        !> Are we in gamma
        logical, intent(in) :: is_Gamma_point

        ! Clean previous
        OMP_OFFLOAD target exit data map(always, delete: polarizability) if(allocated(polarizability))
        if (allocated(polarizability)) deallocate(polarizability)
        if (allocated(nocoulomb_dielectric_wing1)) deallocate( nocoulomb_dielectric_wing1)
        if (allocated(nocoulomb_dielectric_wing2)) deallocate( nocoulomb_dielectric_wing2)

        ! q-dependent polarizability function
        allocate(polarizability(mbsiz,mbsiz,iomstart:iomend), source=zzero)
        OMP_OFFLOAD target enter data map(always, to: polarizability)

        ! For Gamma also the wings
        if (is_Gamma_point) then
          allocate(nocoulomb_dielectric_wing1(mbsiz,3,iomstart:iomend), source=zzero)
          allocate(nocoulomb_dielectric_wing2(mbsiz,3,iomstart:iomend), source=zzero)
        end if

    end subroutine init_polarizability

    subroutine delete_polarizability()

        if (allocated(polarizability)) then
          OMP_OFFLOAD target exit data map(always, delete: polarizability)
          deallocate(polarizability)
        end if

        if (allocated(nocoulomb_dielectric_wing1)) then
          OMP_OFFLOAD target exit data map(always, delete:  nocoulomb_dielectric_wing1)
          deallocate(nocoulomb_dielectric_wing1)
        end if

        if (allocated(nocoulomb_dielectric_wing2)) then
          OMP_OFFLOAD target exit data map(always, delete:  nocoulomb_dielectric_wing2)
          deallocate(nocoulomb_dielectric_wing2)
        end if

    end subroutine delete_polarizability


    subroutine write_polarizability_to_file( iq, is_Gamma_point, file_format, irreducible)

      integer(i32), intent(in)  :: iq 
      logical, intent(in)       :: is_Gamma_point
      character(len=*), intent(in) :: file_format
      logical, intent(in), optional :: irreducible
    
      character(len=max_string_length) :: file_name

      logical :: irreducible_local 

      OMP_OFFLOAD target update from(polarizability)

      if (present(irreducible)) then
        irreducible_local = irreducible
      else
        irreducible_local = .false.
      end if

      if (irreducible_local) then 
        call build_file_name( file_name_polarizability_irreducible, iq, file_name )
      else
        call build_file_name( file_name_polarizability, iq, file_name )
      end if

      call write_to_file( file_name, polarizability, file_format )

      if( is_Gamma_point ) then
          call build_file_name( file_name_nocoulomb_dielectric_wing1, file_name )
          call write_to_file( file_name,  nocoulomb_dielectric_wing1, file_format )
          call build_file_name( file_name_nocoulomb_dielectric_wing2, file_name )
          call write_to_file( file_name,  nocoulomb_dielectric_wing2, file_format )
      end if
      
  end subroutine write_polarizability_to_file


  subroutine read_polarizability_from_file( iq, is_Gamma_point, file_format, irreducible )

      integer(i32), intent(in)  :: iq 
      logical, intent(in)       :: is_Gamma_point
      character(len=*), intent(in) :: file_format
      logical, intent(in), optional :: irreducible
    
      character(len=max_string_length) :: file_name

      logical :: irreducible_local 

      if (present(irreducible)) then
        irreducible_local = irreducible
      else
        irreducible_local = .false.
      end if

      if (irreducible_local) then
        call build_file_name( file_name_polarizability_irreducible, iq, file_name )
      else 
        call build_file_name( file_name_polarizability, iq, file_name )
      end if

      call read_from_file( file_name, polarizability, file_format )

      OMP_OFFLOAD target update to(polarizability)

      if( is_Gamma_point ) then
          call build_file_name( file_name_nocoulomb_dielectric_wing1, file_name )
          call read_from_file( file_name,  nocoulomb_dielectric_wing1, file_format )
          call build_file_name( file_name_nocoulomb_dielectric_wing2, file_name )
          call read_from_file( file_name,  nocoulomb_dielectric_wing2, file_format )
      end if

  end subroutine

  !> Computes the polarizability for a given q-point. For q=0 it computes only
  !> the body, and the dielectric function wings without the Coulomb potential.
  !> The head is ommited, as it is 0, and does not contain
  !> expansion coefficients.
  subroutine compute_polarizability_at_q(iq, iomstart, iomend, Gamma)

      use constants, only: zone
      use modinput, only: input
      use mod_head_and_wings, only: calcwings
      use mod_dielectric_function, only: pmatvv, pmatcv
      use mod_APW_LO, only: apwordmax
      use mod_atoms, only: natmtot
      use mod_muffin_tin, only: lmmaxapw
      use mod_core_states,   only: ncg
      use modgw, only: kqset, fnm, fnm_tet, fnm_sum, mblksiz, b2mb, gkqset, msize
      use mod_product_basis, only: matsiz, minmmat
      use mod_bands, only: nstdf, nomax, numin, eveckpalm, eveckalm, eveck, eveckp
      use mod_eigenvalue_occupancy, only: nstfv
      use mod_eigensystem, only: nmatmax 
      use iso_c_binding,         only: c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
      use mod_device_offload,    only: device_world
      use mod_pointer_remapping, only: remap_fortran_pointer
      use device_linalg_common_interface, only: zgemm_gpu
      use mod_expand_products,   only: expand_products_generic
      use m_memory_device,       only: allocate_device_memory, deallocate_device_memory, &
                                       bytes_double_complex, bytes_int, get_device_pointer


      ! input/output
      integer(i32), intent(in) :: iq
      integer(i32), intent(in) :: iomstart, iomend
      logical, intent(in) :: Gamma
      ! local
      integer(i32) :: ie1, ie2, ibasis
      integer(i32) :: iom
      integer(i32) :: ik, jk, ispn
      integer(i32) :: im, iop, jop
      integer(i32) :: ndim, mdim, nmdim
      integer(i32) :: nblk, iblk, mstart, mend
      real(dp)    :: tstart, tend
      real(dp)    :: wto, wlo
      type(c_ptr) :: minm_cptr
      complex(dp), pointer, contiguous :: minm(:,:,:)
      complex(dp), allocatable :: evecfv(:,:)
      integer(i32) :: my_device
      logical :: tetrahedron_method

      !=============================
      ! Initialization
      !=============================

      ! Get device id
      my_device = device_world%get_device()

      ! total number of states including the core ones
      if (input%gw%coreflag=='all') then
          ndim = nomax+ncg
      else
          ndim = nomax
      end if
      mdim = nstdf-numin+1

      ! block size
      if (mblksiz >= mdim) then
        nblk = 1
      else
        nblk = mdim / mblksiz
        if ( mod(mdim,mblksiz) /= 0 ) nblk = nblk+1
      end if

      ! arrays to store products of KS eigenvectors with the matching coefficients
      allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
      allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
      allocate(eveck(nmatmax,nstfv))
      allocate(eveckp(nmatmax,nstfv))

      OMP_OFFLOAD target enter data map(alloc: eveckalm, eveckpalm, eveck, eveckp)

      !==================================================
      ! Calculate the q-dependent BZ integration weights
      !==================================================
      select case (trim(input%gw%qdepw))
        case('sum')
          tetrahedron_method = .false.
          allocate( fnm_sum(ndim, numin:nstdf, iomstart:iomend), source=zzero )
        case('tet')
          tetrahedron_method = .true.
          allocate( fnm_tet(ndim, numin:nstdf, iomstart:iomend, kqset%nkpt), source=zzero )
          call qdepwtet(iq, iomstart, iomend, ndim)
        case default
          call terminate( "Error(compute_polarizability_at_q): Unknown qdepw method!" )
      end select
  
      !=================
      ! BZ integration
      !=================
      do ik = 1, kqset%nkpt
        if( tetrahedron_method ) then 
          fnm(1:ndim, numin:nstdf, iomstart:iomend) => fnm_tet(:, :, :, ik)
        else
          call qdepwsum(iq, ik, iomstart, iomend, ndim)
          fnm(1:ndim, numin:nstdf, iomstart:iomend) => fnm_sum
        end if

        ! k-q point
        jk = kqset%kqid(ik, iq)

        if (Gamma) call getpmatkgw(ik)

        ! get KS eigenvectors
        allocate(evecfv(nmatmax,nstfv))
        call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
        eveckp = conjg(evecfv)
        call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
        eveck = evecfv
        deallocate(evecfv)

        ! compute products \sum_G C_{k}n * A_{lm}
        call expand_evec(ik,'t')
        call expand_evec(jk,'c')

        OMP_OFFLOAD target update to(eveck, eveckp, eveckalm, eveckpalm)

        !=================================================
        ! Loop over m-blocks in M^i_{nm}(\vec{k},\vec{q})
        !=================================================
        do iblk = 1, nblk

            mstart = numin + (iblk-1)*mblksiz
            mend   = min(nstdf, mstart+mblksiz-1)
            nmdim  = ndim * (mend-mstart+1)

            allocate(minmmat(matsiz,ndim,mstart:mend))
            OMP_OFFLOAD target enter data map(alloc: minmmat)
            msize = sizeof(minmmat)*b2mb

            ! Compute the naked M^i_{nm}+M^i_{cm}, i.e. the ones not multiplied by the barce Coulomb potential
            call expand_products_generic(ik, iq, 1, nomax, 1, ndim-nomax, mstart, mend, 1, 0, minmmat, .false.)

            if (Gamma) then
                ! Wings of the dielectric matrix
                ! Note that the multiplication by the Coulomb potential is missing in this case
                call calcwings(ik, iq, iomstart, iomend, ndim, mstart, mend, matsiz, minmmat,  nocoulomb_dielectric_wing1,  nocoulomb_dielectric_wing2)
            end if

            ! Body
            call allocate_device_memory(minm_cptr, matsiz*nmdim*bytes_double_complex, my_device)
            call c_f_pointer(minm_cptr, minm, int([matsiz,ndim,(mend-mstart+1)],kind=c_size_t))
            ! Remapping the pointer boundaries from Fortran default
            ! TODO(mrm): When supported use lower for c_f_pointer introduced in Fortran 2023
            call remap_fortran_pointer(minm, int([1, 1, mstart], kind=i32), int([matsiz, ndim, mend], kind=i32))

            do iom = iomstart, iomend
#if defined(FLANG_OPENMP_SLICE_MAP_BUG_WORKAROUND)
                OMP_OFFLOAD target data map(to: fnm)
#else
                OMP_OFFLOAD target data map(to: fnm(:,mstart:mend,iom))
#endif
                OMP_OFFLOAD target has_device_addr(minm)
                !$omp teams distribute parallel do collapse(3) default(none) private(ie1,ie2,ibasis) &
                !$omp shared(mstart,mend,ndim,matsiz,minm,fnm,minmmat,iom)
                do ie2 = mstart, mend
                    do ie1 = 1, ndim
                        do ibasis = 1, matsiz
                            minm(ibasis,ie1,ie2) = fnm(ie1,ie2,iom) * &
                                                   minmmat(ibasis,ie1,ie2)
                        end do
                    end do ! ie1
                end do ! ie2
                !$omp end teams distribute parallel do
                OMP_OFFLOAD end target
                OMP_OFFLOAD end target data
                call zgemm_gpu( 'n', 'c', matsiz, matsiz, nmdim, &
                            zone, minm_cptr, matsiz, get_device_pointer(minmmat,my_device), matsiz, &
                            zone, get_device_pointer(polarizability(1,1,iom),my_device), matsiz, device_world)
                call device_world%synchronize()
            end do ! iom
            nullify(minm)
            call deallocate_device_memory(minm_cptr,my_device)
            OMP_OFFLOAD target exit data map(delete: minmmat)
            deallocate(minmmat)

        end do ! iblk

    end do ! ik

    !-------------------
    ! Clear memory
    !-------------------
    if (Gamma) then
        ! deallocate the momentum matrix elements
        deallocate(pmatvv)
        if (input%gw%coreflag=='all') deallocate(pmatcv)
    end if

    OMP_OFFLOAD target update from(polarizability)

    OMP_OFFLOAD target exit data map(delete: eveck, eveckp, eveckalm, eveckpalm)

    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    if( allocated(fnm_sum) ) deallocate( fnm_sum )
    if( allocated(fnm_tet) ) deallocate( fnm_tet )
    if( associated(fnm) ) nullify( fnm )

  end subroutine compute_polarizability_at_q

  !> Computes the dielectric matrix from the polarizability
  subroutine from_polarizability_to_epsilon(iq, Gamma, iomstart, iomend)

      use constants, only: zone, zzero, zi
      use modmpi,    only: rank
      use modxs,      only : symt2
      use modinput, only: input
      use mod_head_and_wings, only: calchead
      use mod_dielectric_function, only: pmatvv, pmatcv, epsh, epsw1, epsw2, epsilon
      use mod_core_states,   only: ncg
      use modgw, only: kqset, fnm, fnm_sum, fnm_tet, freq
      use mod_product_basis, only: matsiz, mbsiz
      use mod_bands, only: nomax, numin, nstdf
      use mod_coulomb_potential, only: barc
      use iso_c_binding,         only: c_ptr, c_loc, c_f_pointer, c_sizeof, c_size_t
      use mod_device_offload,    only: device_world
      use mod_pointer_remapping, only: remap_fortran_pointer
      use device_linalg_common_interface, only: zgemm_batched_gpu, zgemm_gpu
      use m_memory_device,       only: allocate_device_memory, deallocate_device_memory, &
                                       bytes_double_complex, bytes_int, get_device_pointer

      ! input/output
      integer(i32), intent(in) :: iq
      logical, intent(in)      :: Gamma
      integer(i32), intent(in) :: iomstart, iomend
      ! local
      integer(i32) :: ie1, ie2
      integer(i32) :: iom
      integer(i32) :: ik
      integer(i32) :: im, iop, jop
      integer(i32) :: ndim
      integer(i32) :: nblk, iblk, mstart, mend
      real(dp)    :: wto, wlo
      complex(dp) :: head(3,3), f, w
      complex(dp), allocatable :: evecfv(:,:)
      type(c_ptr) :: temp_vcpol_cptr
      integer(i32) :: my_device
      logical :: tetrahedron_method

      my_device = device_world%get_device()

      ! total number of states including the core ones
      if (input%gw%coreflag=='all') then
          ndim = nomax+ncg
      else
          ndim = nomax
      end if

      ! For Gamma point compute the head of the dielectric matrix
      if (Gamma) then

        ! Calculate the q-dependent BZ integration weights
        select case (trim(input%gw%qdepw))
          case('sum')
            tetrahedron_method = .false.
            allocate( fnm_sum(ndim, numin:nstdf, iomstart:iomend), source=zzero )
          case('tet')
            tetrahedron_method = .true.
            allocate( fnm_tet(ndim, numin:nstdf, iomstart:iomend, kqset%nkpt), source=zzero )
            call qdepwtet(iq, iomstart, iomend, ndim)
          case default
            call terminate( "Error(compute_polarizability_at_q): Unknown qdepw method!" )
        end select

        do ik = 1, kqset%nkpt
            if( tetrahedron_method ) then 
              fnm(1:ndim, numin:nstdf, iomstart:iomend) => fnm_tet(:, :, :, ik)
            else
              call qdepwsum(iq, ik, iomstart, iomend, ndim)
              fnm(1:ndim, numin:nstdf, iomstart:iomend) => fnm_sum
            end if
            ! Read the momentum matrix elements
            call getpmatkgw(ik)
            ! Compute the head of the dielectric function
            call calchead(ik, numin, nstdf, iomstart, iomend, ndim, epsh)

        end do

        ! Clear memory
        deallocate(pmatvv)
        if (input%gw%coreflag=='all') deallocate(pmatcv)
        if( allocated(fnm_sum) ) deallocate( fnm_sum )
        if( allocated(fnm_tet) ) deallocate( fnm_tet )
        if( associated(fnm) ) nullify( fnm )


        ! Compute the wings of the dielectric matrix
        OMP_OFFLOAD target data map(to: nocoulomb_dielectric_wing1, nocoulomb_dielectric_wing2)
        call zgemm_batched_gpu('c', 'n', mbsiz, 3, matsiz, zone, &
                    get_device_pointer(barc,my_device), matsiz, &
                    get_device_pointer(nocoulomb_dielectric_wing1,my_device), matsiz,  &
                    zzero, get_device_pointer(epsw1(1,1,iomstart),my_device), mbsiz, &
                    iomend-iomstart+1, device_world, stridea=0_i32)
        call device_world%synchronize()
        call zgemm_batched_gpu('t', 'n', mbsiz, 3, matsiz, zone, &
                    get_device_pointer(barc,my_device), matsiz, &
                    get_device_pointer(nocoulomb_dielectric_wing2,my_device), matsiz,  &
                    zzero, get_device_pointer(epsw2(1,1,iomstart),my_device), mbsiz, &
                    iomend-iomstart+1, device_world, stridea=0_i32)
        call device_world%synchronize()
        OMP_OFFLOAD end target data
      
      end if

      ! Allocating the temporary. For device aware compilation it does live in the device
      call allocate_device_memory(temp_vcpol_cptr, matsiz*mbsiz*bytes_double_complex, my_device)
      do iom = iomstart, iomend
          call zgemm_gpu('c', 'n', mbsiz, matsiz, matsiz, zone, &
                        get_device_pointer(barc,my_device), matsiz, &
                        get_device_pointer(polarizability(1,1,iom),my_device), matsiz, &
                        zzero, temp_vcpol_cptr, mbsiz, device_world)
          call device_world%synchronize()
          call zgemm_gpu('n', 'n', mbsiz, mbsiz, matsiz, zone, &
                        temp_vcpol_cptr, mbsiz, &
                        get_device_pointer(barc,my_device), matsiz, &
                        zzero, get_device_pointer(epsilon(1,1,iom),my_device), mbsiz, device_world)
          call device_world%synchronize()
      end do
      call deallocate_device_memory(temp_vcpol_cptr, my_device)

      OMP_OFFLOAD target update from(epsilon)

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

  end subroutine from_polarizability_to_epsilon

end module
