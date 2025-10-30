!> Module handling electronic part of an EPH calculation.
module eph_electrons
  use eph_variables

  use precision, only: dp
  use asserts, only: assert
  use exciting_mpi, only: xmpi_allgatherv
  use modmpi
  use matrix_fourier_interpolation, only: mfi_type
  use block_data_file, only: block_data_file_type

  implicit none
  private

  !> electron energies on electron \({\bf k}\)-grid (relative to Fermi level)
  real(dp), allocatable, public :: eph_el_energy_k(:,:)
  !> Wannier transformation matrices \(U_{mn}({\bf k})\) on electron \({\bf k}\)-grid
  complex(dp), allocatable, public :: eph_el_evec_k(:,:,:)
  !> object for matrix Fourier interpolation on electron \({\bf k}\)-grid
  type(mfi_type), public :: eph_el_mfi
  !> electron Hamiltonian in real space Wannier gauge, \(\mathcal{H}_{mn}({\bf R})\)
  complex(dp), allocatable :: eph_el_HR(:,:,:)
  !> name for binary file to save Hamiltonian in real space Wannier gauge for later access
  character(*), parameter :: eph_el_HR_filename = "EPH_HR.OUT"
  !> name for binary file to save Wannier gauge matrices for later access
  character(*), parameter :: eph_el_Uk_filename = "EPH_Uk.OUT"

  public :: eph_el_free, eph_el_set_energies, eph_el_setup_interpolation, eph_el_interpolate

contains

  !> Free memory from module variables.
  subroutine eph_el_free
    if (allocated(eph_el_energy_k)) deallocate( eph_el_energy_k )
    if (allocated(eph_el_evec_k)) deallocate( eph_el_evec_k )
    if (allocated(eph_el_HR)) deallocate( eph_el_HR )
    call eph_el_mfi%destroy
  end subroutine eph_el_free

  !================================================================================ 
  ! SET CORRECT ELECTRON ENERGIES
  !
  !> Replace KS energies by the ones used in the Wannier calculation according to
  !> the attribute `input` in `<wannier>`.
  !>
  !> Throughout the eph calculation, we will access energies via [[dfpt_eig_geteval(subroutine)]].
  !> We replace the eigenenergies in the respective file by the ones obtained from [[wfhelp_geteval(subroutine)]].
  subroutine eph_el_set_energies
    use dfpt_variables, only : dfpt_kset, fevalk0
    use mod_wannier_variables, only : wf_kset
    use mod_wannier_helper, only : wfhelp_geteval
    use modinput

    integer :: fst, lst, ik_dfpt, ik_wan, isym

    integer, allocatable :: shp(:)
    real(dp), allocatable :: eval_wan(:,:), eval_dfpt(:)

    ! return, if there is nothing to do
    if (.not. eph_use_wannier) return
    if (input%properties%wannier%input == 'gs') return

    ! read energies used for Wannier functions
    call wfhelp_geteval( eval_wan, fst, lst )
    ! read and replace DFPT energies
    shp = fevalk0%get_block_shape()
    allocate( eval_dfpt(shp(1)) )
    do ik_dfpt = 1, dfpt_kset%nkpt
      call fevalk0%read( ik_dfpt, eval_dfpt )
      call findkptinset( dfpt_kset%vkl(:, ik_dfpt), wf_kset, isym, ik_wan )
      eval_dfpt(fst:lst) = eval_wan(:, ik_wan)
      ! shift energies relative to Fermi level and apply scissor
      eval_dfpt = eval_dfpt - eph_efermi
      where( eval_dfpt < 0.0_dp )
        eval_dfpt = eval_dfpt - eph_scissor(1)
      elsewhere
        eval_dfpt = eval_dfpt + eph_scissor(2)
      end where
      call fevalk0%write( ik_dfpt, eval_dfpt )
    end do
    deallocate( eval_wan, eval_dfpt, shp )
  end subroutine eph_el_set_energies
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! SETUP WANNIER INTERPOLATION OF ELECTRONS
  !
  !> Set up the Wannier interpolation of the electron Hamiltonian.
  !>
  !> This includes
  !> 
  !>   * checking for availability of Wannier functions
  !>   * preparation of matrix Fourier interpolation on electron \({\bf k}\)-grid
  !>   * computation of localized Hamiltonian in real space Wannier gauge \(\mathcal{H}_{mn}({\bf R})\)
  !>     and writing to file
  !>   * or reading \(\mathcal{H}_{mn}({\bf R})\) from file, if possible
  subroutine eph_el_setup_interpolation( write_localization )
    use dfpt_variables, only: dfpt_Gset, dfpt_Gkset
    use mod_kpointset, only: generate_Gk_vectors
    !> write spatial localization of \(\mathcal{\bf H}({\bf R})\) to file (default: `.false.`)
    logical, optional, intent(in) :: write_localization

    integer :: ik, ik0, ir, isym, un, stat
    logical :: write_loc
    type(block_data_file_type) :: HR_file

    real(dp), allocatable :: centers(:,:)
    complex(dp), allocatable :: Hk(:,:,:)

    write_loc = .false.
    if (present(write_localization)) write_loc = write_localization

    ! delete existing module variables
    call eph_el_free

    ! try to read Wannier functions and electron energies on support grid
    if (.not. allocated(eph_Gkset_el%ngk)) &
      call generate_Gk_vectors( eph_Gkset_el, eph_kset_el, dfpt_Gset, dfpt_Gkset%gkmax )
    call eph_el_read_wannier( eph_kset_el, eph_Gkset_el, eph_nwf_tot, eph_el_energy_k, eph_el_evec_k, centers )
    call terminate_if_false( eph_nwf_tot > 0, '(eph_el_setup_interpolation) &
      No Wannier functions could be found. Make sure, they have been precomputed.' )

    ! shift energies relative to Fermi level and apply scissor
    eph_el_energy_k = eph_el_energy_k - eph_efermi
    where( eph_el_energy_k < 0.0_dp )
      eph_el_energy_k = eph_el_energy_k - eph_scissor(1)
    elsewhere
      eph_el_energy_k = eph_el_energy_k + eph_scissor(2)
    end where

    ! setup matrix fourier interpolation
    ! (Note: The electron / Wannier k-grid is assumed to be non-reduced!)
    eph_el_mfi = mfi_type( eph_kset_el%ngridk, eph_kset_el%bvec, eph_kset_el%vkl )
    call eph_el_mfi%set_localization_centers( left_centers=centers, right_centers=centers )

    ! set up binary file
    HR_file = block_data_file_type( eph_el_HR_filename, [eph_nwf_tot, eph_nwf_tot], cmplx( 0, 0, dp ) )

    ! allocate real space Wannier Hamiltonian
    allocate( eph_el_HR(eph_nwf_tot, eph_nwf_tot, eph_el_mfi%nr) )

    ! try to read H_W(R) from file
    if (HR_file%exists()) then
      call HR_file%open( mpiglobal )
      do ir = 1, eph_el_mfi%nr
        call HR_file%read( ir, eph_el_HR(:, :, ir) )
      end do
      call HR_file%close( mpiglobal )
    ! compute H_W(R) and write to file
    else
      allocate( Hk(eph_nwf_tot, eph_nwf_tot, eph_el_mfi%np) )
      do ik = 1, eph_el_mfi%np
        call findkptinset( eph_el_mfi%vpl(:, ik), eph_kset_el, isym, ik0 )
        call gen_Hk_wannier( eph_el_energy_k(:, ik0), eph_el_evec_k(:, :, ik0), Hk(:, :, ik) )
      end do
      call eph_el_mfi%transform_p2R( [eph_nwf_tot, eph_nwf_tot], 1, Hk, eph_nwf_tot**2, 1, eph_el_HR, eph_nwf_tot**2, 1 )
      deallocate( Hk )

      call HR_file%open( mpiglobal )
      do ir = 1, eph_el_mfi%nr
        call HR_file%write( ir, eph_el_HR(:, :, ir) )
      end do
      call HR_file%close( mpiglobal )
    end if

    if (allocated(centers)) deallocate( centers )

    ! write spatial localization to file
    if (mpiglobal%rank == 0 .and. write_loc) then
      open( newunit=un, file='eph_el_loc.dat', action='write', form='formatted', iostat=stat )
      call terminate_if_false( stat == 0, '(eph_el_setup_interpolation) &
        Failed to open file `eph_el_loc.dat`.' )
      write( un, '("#",a5,a26,a26)' ) 'iR', '|R|', 'max(|H(R)|)'
      do ir = 1, eph_el_mfi%nr
        write( un, '(i6,2g26.16)' ) ir, eph_el_mfi%rlen(ir), maxval( abs( eph_el_HR(:, :, ir) ) )
      end do
      close( un )
    end if
  end subroutine eph_el_setup_interpolation
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! INTERPOLATE ELECTRONS
  !
  !> Get electron energies \(\epsilon_{n}({\bf k}')\) and eigenvectors \(U_{mn}({\bf k}')\)
  !> for a given set of wave vectors \({\bf k}'\) by Wannier interpolation.
  !>
  !> This is done by Fourier interpolating the Hamiltonian to \({\bf k}'\) by
  !> \[ \mathcal{H}({\bf k}') = \sum_{\bf R} {\rm e}^{{\rm i} {\bf k}'\cdot{\bf R}}\, \mathcal{H}({\bf R}) \]
  !> and diagonalizing \(\mathcal{H}({\bf k}')\) such that
  !> \[ H({\bf k}') = \operatorname{diag}\left(\epsilon_{n}({\bf k}')\right) 
  !>    = U({\bf k}')\, \mathcal{H}({\bf k}')\, U^\dagger({\bf k}') \;. \]
  !> See also [[transform_R2p(subroutine)]]. Electron energies are returned relative to the Fermi level.
  !>
  !> MPI parallelization is over \({\bf k}'\) points, but each process can specify a different band range `irange`.
  subroutine eph_el_interpolate( vkl, evalk, Umnk, &
      irange, mpicomm )
    use exciting_mpi, only: xmpi_allreduce
#ifdef MPI
    use mpi_f08, only: MPI_Send, MPI_Recv, MPI_DOUBLE, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, MPI_COMM
#endif
    use m_linalg, only: zhediag
    !> set of wave vectors \({\bf k}'\) in lattice coordinates
    real(dp), intent(in) :: vkl(:,:)
    !> electron energies \(\epsilon_{n}({\bf k}')\)
    real(dp), allocatable, intent(out) :: evalk(:,:)
    !> electron eigenvectors \(U_{mn}({\bf k}')\)
    complex(dp), allocatable, intent(out) :: Umnk(:,:,:)
    !> range of bands for which to interpolate (default: `[eph_fwf, eph_lwf]`)
    integer, optional, intent(in) :: irange(2)
    !> MPI communicator (default: global MPI communicator)
    type(mpiinfo), optional, intent(inout) :: mpicomm

    integer :: nk, ik, ik1, ik2, nwf, irng(2), irank, jrank
    type(mpiinfo) :: mpi

    integer, allocatable :: irng_list(:,:), ik_range(:,:), ix_range(:,:)
    real(dp), allocatable :: eval(:,:)
    complex(dp), allocatable :: H(:,:,:), evec(:,:)

    irng = [eph_fwf, eph_lwf]
    if (present(irange)) irng = irange

    mpi = mpiglobal
    if (present(mpicomm)) mpi = mpicomm

    nwf = irng(2) - irng(1) + 1
    nk = size( vkl, dim=2 )

    call assert( size( vkl, dim=1 ) == 3, &
      '`vkl` must be a set of vectors of length 3.' )

    ! communicate information on band distribution
    allocate( irng_list(2, mpi%procs) )
    irng_list = 0; irng_list(:, mpi%rank+1) = irng
    call xmpi_allreduce( irng_list, mpi )
    
    ! distribute k among processes
    call patchwork_distribution( nk, 1, mpi%procs, ik_range, ix_range )
    ik1 = ik_range(1, mpi%rank+1)
    ik2 = ik_range(2, mpi%rank+1)
    allocate( H(eph_nwf_tot, eph_nwf_tot, ik1:ik2) )

    ! allocate output
    if (allocated(evalk)) deallocate( evalk )
    allocate( evalk(irng(1):irng(2), nk) )
    if (allocated(Umnk)) deallocate( Umnk )
    allocate( Umnk(irng(1):irng(2), eph_nwf_tot, nk) )

    ! interpolate Hamiltonian
    call eph_el_mfi%transform_R2p( [eph_nwf_tot, eph_nwf_tot], 1, eph_el_HR, eph_nwf_tot**2, 1, H, eph_nwf_tot**2, 1, vkl(:, ik1:ik2), &
      minimal_distances=.true. )

    ! diagonalize Hamiltonian
    allocate( eval(eph_nwf_tot, ik1:ik2), evec(eph_nwf_tot, eph_nwf_tot) )
    do ik = ik1, ik2
      call zhediag( H(:, :, ik), eval(:, ik), evec )
      ! U(k) is the h.c. of the eigenvectors
      H(:, :, ik) = conjg( transpose( evec ) )
      evalk(:, ik) = eval(irng(1):irng(2), ik)
      Umnk(:, :, ik) = H(irng(1):irng(2), :, ik)
    end do

    ! distribute results
#ifdef MPI
    do jrank = 0, mpi%procs-1
      if (mpi%rank == jrank) then
        do irank = 0, mpi%procs-1
          if (irank == jrank) cycle
          do ik = ik_range(1, irank+1), ik_range(2, irank+1)
            call MPI_Recv( evalk(:, ik), nwf, MPI_DOUBLE, irank, irank*nk+ik, MPI_COMM(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
            call MPI_Recv( Umnk(:, :, ik), nwf*eph_nwf_tot, MPI_DOUBLE_COMPLEX, irank, (mpi%procs+irank)*nk+ik, MPI_COMM(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
          end do
        end do
      else
        do ik = ik1, ik2
          evec = H(irng_list(1, jrank+1):irng_list(2, jrank+1), :, ik)
          call MPI_Send( evalk(irng_list(1, jrank+1), ik), size( evec, dim=1 ), MPI_DOUBLE, jrank, mpi%rank*nk+ik, MPI_COMM(mpi%comm), mpi%ierr )
          call MPI_Send( evec, size( evec ), MPI_DOUBLE_COMPLEX, jrank, (mpi%procs+mpi%rank)*nk+ik, MPI_COMM(mpi%comm), mpi%ierr )
        end do
      end if
    end do 
#endif
    deallocate( H, eval, evec, irng_list, ik_range, ix_range )
  end subroutine eph_el_interpolate
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! AUXILIARY PROCEDURES
  !
  !> Try to receive precomputed Wannier functions.
  !> This includes the electron energies used for generating the Wannier functions
  !> (might be KS, hybrid or GW energies) and unitary transformations \({\bf U}({\bf k})\)
  !> defining the Wannier functions.
  !>
  !> The unitary transformations are rotated to match the phase of the eigenvectors from the
  !> DFPT calculation which are used to compute matrix elements.
  subroutine eph_el_read_wannier( kset, Gkset, nwf, evalk, eveck, centers )
    use constants, only: zzero, zone, twopi
    use dfpt_variables, only : dfpt_kset, dfpt_Gkset, feveck0
    use dfpt_eigensystem, only : dfpt_eig_getevec
    use mod_wannier_variables, only: wf_fst, wf_lst, wf_nst, wf_nwf, wf_kset, wf_transform, wf_centers
    use mod_wannier_filehandling, only : wffile_readtransform
    use mod_wannier_helper, only : wfhelp_geteval, wfhelp_getevec
    use mod_kpointset, only: k_set, Gk_set
    use mod_eigensystem, only : nmatmax_ptr
    use mod_eigenvalue_occupancy, only : nstfv
    use mod_spin, only : nspinor
    use mod_APW_LO, only : nlotot
    use m_linalg, only : zlsp
    use xlapack, only : svd_divide_conquer
    !> set of \({\bf k}\)-vectors for which electron energies \(\epsilon_{n{\bf k}}\) 
    !> and Wannier matrices \(U_{mn}({\bf k})\) should be read
    type(k_set), intent(in) :: kset
    !> set of \({\bf G+k}\)-vectors
    type(Gk_set), intent(in) :: Gkset
    !> number of Wannier functions (returns `0`, if no Wannier functions were found)
    integer, intent(out) :: nwf
    !> electron eigenenergies \(\epsilon_{n{\bf k}}\) for all wannierized bands
    real(dp), allocatable, intent(out) :: evalk(:,:)
    !> Wannier matrices \(U_{mn}({\bf k})\) for all wannierized bands and all Wannier functions
    complex(dp), allocatable, intent(out) :: eveck(:,:,:)
    !> localization centers of all Wannier functions lattice coordinates
    real(dp), allocatable, intent(out) :: centers(:,:)

    integer :: fst, lst, ik, ik1, ik2, ikw, isym, nmat
    integer, target :: nmatmax
    logical :: success
    type(block_data_file_type) :: Uk_file

    real(dp), allocatable :: eval(:,:), sval(:)
    complex(dp), allocatable :: evec_wan(:,:,:), evec_dfpt(:,:), rot(:,:), lsvec(:,:), rsvec(:,:)

    nwf = 0

    ! try to read electron energies
    call wfhelp_geteval( eval, fst, lst )

    ! try to read Wannier functions
    call wffile_readtransform( success )
    if (.not. success) return
    
    ! set number of Wannier functions
    nwf = wf_nwf

    ! read energies
    if (allocated(evalk)) deallocate( evalk )
    allocate( evalk(wf_fst:wf_lst, kset%nkpt) )
    do ik = 1, kset%nkpt
      call findkptinset( kset%vkl(:, ik), wf_kset, isym, ikw )
      call terminate_if_false( isym == 1, '(eph_el_read_wannier) &
        Requested k-point not found in Wannier k-point set.' )
      evalk(:, ik) = eval(wf_fst:wf_lst, ikw)
    end do

    ! get Wannier gauge matrices U(k)
    if (allocated(eveck)) deallocate( eveck )
    allocate( eveck(wf_fst:wf_lst, wf_nwf, kset%nkpt) )
    Uk_file = block_data_file_type( eph_el_Uk_filename, [wf_nst, wf_nwf], cmplx( 0, 0, dp ) )
    ik1 = firstofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
    ik2 = lastofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
    ! try to read rotated U(k) from file
    if (Uk_file%exists()) then
      call Uk_file%open( mpiglobal )
      do ik = ik1, ik2
        call Uk_file%read( ik, eveck(:, :, ik) )
      end do
      call Uk_file%close( mpiglobal )
    ! read and rotate original U(k) and write to file
    else
      nmatmax = Gkset%ngkmax + nlotot
      nmatmax_ptr => nmatmax
      allocate( evec_wan(nmatmax_ptr, nstfv, nspinor) )

      call Uk_file%open( mpiglobal )
      do ik = ik1, ik2
        nmat = Gkset%ngk(1, ik) + nlotot
        allocate( rot(nmat, wf_nst), sval(wf_nst), lsvec(nmat, wf_nst), rsvec(wf_nst, wf_nst) )
        call findkptinset( kset%vkl(:, ik), wf_kset, isym, ikw )
        call wfhelp_getevec( ikw, evec_wan )
        call dfpt_eig_getevec( kset%vkl(:, ik), Gkset%vgkl(:, :, 1, ik), feveck0, dfpt_kset, dfpt_Gkset, [1, nmat], evec_dfpt )
        ! find matrix that transforms between DFPT and Wannier eigenvectors ...
        call zlsp( evec_dfpt(:nmat, :), evec_wan(:nmat, wf_fst:wf_lst, 1), rot )
        ! ... and make it unitary
        call svd_divide_conquer( rot, sval, U=lsvec, V_H=rsvec )
        call zgemm( 'n', 'n', nmat, wf_nst, wf_nst, zone, lsvec, nmat, rsvec, wf_nst, zzero, rot, nmat )
        call zgemm( 'n', 'n', wf_nst, wf_nwf, wf_nst, zone, &
          rot(wf_fst, 1), nmat, &
          wf_transform(:, :, ikw), wf_nst, zzero, &
          eveck(:, :, ik), wf_nst )
        deallocate( rot, sval, lsvec, rsvec )
        call Uk_file%write( ik, eveck(:, :, ik) )
      end do
      call Uk_file%close( mpiglobal )

      deallocate( evec_wan )
      if (allocated(eval)) deallocate( eval )
      if (allocated(evec_dfpt)) deallocate( evec_dfpt )
    end if
    call xmpi_allgatherv( mpiglobal, eveck, wf_nst * wf_nwf * (ik2 - ik1 + 1) )

    ! read localization centers (transformed to lattice coordinates)
    if (allocated(centers)) deallocate( centers )
    allocate( centers(3, wf_nwf) )
    call dgemm( 't', 'n', 3, wf_nwf, 3, 1.0_dp/twopi, wf_kset%bvec, 3, wf_centers, 3, 0.0_dp, centers, 3 )
  end subroutine eph_el_read_wannier

  !> Generate the Hamiltonian at \({\bf k}\) in Wannier gauge, \(\mathcal{H}_{mn}({\bf k})\).
  subroutine gen_Hk_wannier( evalk, Umnk, Hk )
    use constants, only: zzero, zone
    use math_utils, only: is_square
    !> electron energies at \({\bf k}\) (relative to Fermi energy)
    real(dp), intent(in) :: evalk(:)
    !> Wannier transformation matrix \(U_{mn}({\bf k})\)
    complex(dp), intent(in) :: Umnk(:,:)
    !> Hamiltonian in Wannier gauge \(\mathcal{H}_{mn}({\bf k})\)
    complex(dp), intent(out) :: Hk(:,:)

    integer :: nst, nwf, ist

    complex(dp), allocatable :: auxmat(:,:)

    nst = size( Umnk, dim=1 )
    nwf = size( Umnk, dim=2 )

    call assert( is_square( Hk ), &
      '`Hk` must be square.' )
    call assert( nst == size( evalk ), &
      'Number of electronic states in `Umnk` and `evalk` do not match.' )
    call assert( nwf == size( Hk, dim=1 ), &
      'Number of Wannier functions in `Umnk` and `Hk` do not match.' )

    allocate( auxmat(nst, nwf) )
    do ist = 1, nst
      auxmat(ist, :) = evalk(ist) * Umnk(ist, :)
    end do
    call zgemm( 'c', 'n', nwf, nwf, nst, zone, Umnk, nst, auxmat, nst, zzero, Hk, nwf )
    deallocate( auxmat )
  end subroutine gen_Hk_wannier
  !-------------------------------------------------------------------------------- 
end module eph_electrons
