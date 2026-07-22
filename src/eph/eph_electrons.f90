!> Module handling electronic part of an EPH calculation.
module eph_electrons
  use eph_variables

  use precision, only: dp
#include "asserts.fpp"
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
  !> apply minimal distance interpolation for electrons
  logical, public :: eph_el_mindist = .true.
  !> electron Hamiltonian in real space Wannier gauge, \(\mathcal{H}_{mn}({\bf R})\)
  complex(dp), allocatable :: eph_el_HR(:,:,:)
  !> name for binary file to save Hamiltonian in real space Wannier gauge for later access
  character(*), parameter :: eph_el_HR_filename = "EPH_HR.OUT"
  !> name for binary file to save Wannier gauge matrices for later access
  character(*), parameter :: eph_el_Uk_filename = "EPH_Uk.OUT"

  public :: eph_el_free, eph_el_set_wannier_eigensystem, eph_el_setup_interpolation, eph_el_interpolate, eph_el_gen_Hk_wannier, eph_el_set_default_frequency_grid, eph_el_fermi_and_scissor

contains

  !> Free memory from module variables.
  subroutine eph_el_free
    if (allocated(eph_el_energy_k)) deallocate( eph_el_energy_k )
    if (allocated(eph_el_evec_k)) deallocate( eph_el_evec_k )
    if (allocated(eph_el_HR)) deallocate( eph_el_HR )
    call eph_el_mfi%destroy
  end subroutine eph_el_free

  !================================================================================ 
  ! SET CORRECT ELECTRON ENERGIES AND EIGENVECTORS
  !
  !> Replace KS energies and eigenvectors by the ones used in the Wannier calculation according to
  !> the attribute `input` in `<wannier>`.
  !>
  !> Throughout the eph calculation, we will access energies and eigenvectors via [[dfpt_eig_geteval(subroutine)]]
  !> and [[dfpt_eig_getevec(subroutine)]], respectively.
  !> We replace the eigenenergies and eigenvectors in the respective file by the ones obtained from 
  !> [[wfhelp_geteval(subroutine)]] and [[wfhelp_getevec(subroutine)]].
  subroutine eph_el_set_wannier_eigensystem
    use dfpt_variables, only: dfpt_kset, dfpt_Gkset, fevalk0, feveck0
    use mod_wannier_variables, only: wf_kset
    use mod_wannier_helper, only: wfhelp_geteval, wfhelp_getevec
    use mod_eigensystem, only: nmatmax_ptr
    use mod_eigenvalue_occupancy, only: nstsv
    use mod_APW_LO, only: nlotot
    use mod_spin, only: nspinor
    use sorting, only: sort_index_1d
    use modinput

    integer :: fst, lst, ik_dfpt, ik_wan, ik1, ik2, isym, nmat, ist
    integer, target :: nmatmax

    integer, allocatable :: sort(:), shp(:)
    real(dp), allocatable :: eval_wan(:,:), eval_dfpt(:)
    complex(dp), allocatable :: evec_wan(:,:,:), evec_dfpt(:,:)

    ! return, if there is nothing to do
    if (.not. eph_use_wannier) return
    if (input%properties%wannier%input == 'gs') return

    ! read energies used for Wannier functions
    call wfhelp_geteval( eval_wan, fst, lst )

    shp = fevalk0%get_block_shape()
    allocate( eval_dfpt(shp(1)) )
    shp = feveck0%get_block_shape()
    allocate( evec_dfpt(shp(1), shp(2)) )
    nmatmax = dfpt_Gkset%ngkmax + nlotot
    nmatmax_ptr => nmatmax
    allocate( evec_wan(nmatmax_ptr, nstsv, nspinor) )

    ! (Resort energies. Might not be sorted in case of GW.)
    sort = [(ist, ist=1, nstsv)]
    ik1 = firstofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
    ik2 = lastofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
    do ik_dfpt = ik1, ik2
      call findkptinset( dfpt_kset%vkl(:, ik_dfpt), wf_kset, isym, ik_wan )
      nmat = dfpt_Gkset%ngk(1, ik_dfpt) + nlotot
      sort(fst:lst) = sort_index_1d( size(eval_wan, dim=1), eval_wan(:, ik_wan) ) + fst - 1
      ! replace DFPT energies
      call fevalk0%read( ik_dfpt, eval_dfpt )
      eval_dfpt(fst:lst) = eval_wan(sort(fst:lst), ik_wan)
      call fevalk0%write( ik_dfpt, eval_dfpt )
      ! replace DFPT eigenvectors
      call wfhelp_getevec( ik_wan, evec_wan )
      call feveck0%read( ik_dfpt, evec_dfpt )
      evec_dfpt(1:nmat, 1:nstsv) = evec_wan(1:nmat, sort, 1)
      call feveck0%write( ik_dfpt, evec_dfpt )
    end do
    deallocate( eval_wan, eval_dfpt, evec_wan, evec_dfpt, shp )
  end subroutine eph_el_set_wannier_eigensystem
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
  !>   * building an index map between Wannier interpolated bands and original input bands
  subroutine eph_el_setup_interpolation( write_localization )
    use dfpt_variables, only: dfpt_Gset, dfpt_Gkset
    use mod_kpointset, only: generate_Gk_vectors
    use modinput
    !> write spatial localization of \(\mathcal{\bf H}({\bf R})\) to file (default: `.false.`)
    logical, optional, intent(in) :: write_localization

    integer :: ik, ik0, ir, isym, ist, jst, un, stat
    logical :: write_loc
    type(block_data_file_type) :: HR_file

    integer, allocatable :: map(:,:)
    real(dp), allocatable :: centers(:,:), eval(:,:)
    complex(dp), allocatable :: Hk(:,:,:), evec(:,:,:)

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
    call eph_el_fermi_and_scissor( eph_el_energy_k, size(eph_el_energy_k), eph_efermi, eph_scissor )

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
        call eph_el_gen_Hk_wannier( eph_el_energy_k(:, ik0), eph_el_evec_k(:, :, ik0), Hk(:, :, ik) )
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

    ! build index map
    ! This assumes a unique and contiguous Wannierization, i.e., the Wannierized part of the
    ! band structure must not contain any gaps and each band must be Wannierized exactly once.
    ! (No two groups with overlapping inner windows.)
    allocate( eph_wf_band_map(eph_nwf_tot, 0:eph_kset_el%nkpt), source=0 )
    call eph_el_interpolate( eph_kset_el%vkl(:, 1:eph_kset_el%nkpt), eval, evec, irange=[1, eph_nwf_tot] )
    do ik = 1, eph_kset_el%nkpt
      jst = eph_fst - 1 + maxloc( [(count( abs( eph_el_energy_k(ist:ist+eph_nwf_tot-1, ik) - eval(:, ik) ) < eph_el_degtol ), ist=eph_fst, eph_lst-eph_nwf_tot+1)], dim=1 )
      eval(:, ik) = eval(:, ik) - eph_el_energy_k(jst:jst+eph_nwf_tot-1, ik) 
      where (abs(eval(:, ik)) < eph_el_degtol)
        eph_wf_band_map(:, ik) = [(ist, ist=jst, jst+eph_nwf_tot-1)]
      end where
    end do
    eph_wf_band_map(:, 0) = sum( eph_wf_band_map, dim=2 ) / eph_kset_el%nkpt
    where (eph_wf_band_map(:, 0) /= maxval(eph_wf_band_map, dim=2)) 
      eph_wf_band_map(:, 0) = 0
    end where
    call terminate_if_false( any( eph_wf_band_map(:, 0) > 0 ), '(eph_el_setup_interpolation) &
      No band seems to be described exactly by Wannierization. Index alignment failed.' )
    if (mpiglobal%rank == 0) then
      write( *, '("Info: Range of Wannier functions that exactly describe original bands is ",i3," to ",i3,".")' ) &
        minloc( pack( eph_wf_band_map(:, 0), eph_wf_band_map(:, 0) > 0 ) ), &
        maxloc( pack( eph_wf_band_map(:, 0), eph_wf_band_map(:, 0) > 0 ) )
      write( *, '("      Range of Wannier functions that are included in EPH calculation is ",i3," to ",i3,".")' ) &
        eph_fwf, eph_lwf
    end if
    ! original index of first and last state that is covered by Wannier functions
    eph_fst_span = findloc( eph_wf_band_map(:, 0) > 0, .true., dim=1 )
    if (eph_fst_span > 0) eph_fst_span = eph_wf_band_map(eph_fwf+eph_fst_span-1, 0) - eph_fst_span + 1
    eph_lst_span = eph_fst_span + eph_nwf - 1
    if (mpiglobal%rank == 0) then
      write( *, '("      Range of original bands that are described by Wannier functions is ",i3," to ",i3,".")' ) &
        eph_fst_span, eph_lst_span
    end if

    ! write spatial localization to file
    call barrier( mpicom=mpiglobal )
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

    integer :: nk, ik, ik1, ik2, nwf, irng(2), irank
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

    CALL_ASSERT( size( vkl, dim=1 ) == 3,  '`vkl` must be a set of vectors of length 3.' )

    ! communicate information on band distribution
    allocate( irng_list(3, 0:mpi%procs-1) )
    irng_list = 0; irng_list(:, mpi%rank) = [irng, irng(2)-irng(1)+1]
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
    if (ik2 >= ik1) &
      call eph_el_mfi%transform_R2p( [eph_nwf_tot, eph_nwf_tot], 1, eph_el_HR, eph_nwf_tot**2, 1, H, eph_nwf_tot**2, 1, vkl(:, ik1:ik2), &
        minimal_distances=eph_el_mindist )

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
    do irank = 0, mpi%procs-1
      if (irank == mpi%rank) cycle
      do ik = ik1, ik2
        evec = H(irng_list(1, irank):irng_list(2, irank), :, ik)
        call MPI_Send( eval(irng_list(1, irank), ik), irng_list(3, irank), MPI_DOUBLE, irank, mpi%rank*nk+ik, MPI_COMM(mpi%comm), mpi%ierr )
        call MPI_Send( evec, size( evec ), MPI_DOUBLE_COMPLEX, irank, (mpi%procs+mpi%rank)*nk+ik, MPI_COMM(mpi%comm), mpi%ierr )
      end do
    end do
    do irank = 0, mpi%procs-1
      if (irank == mpi%rank) cycle
      do ik = ik_range(1, irank+1), ik_range(2, irank+1)
        call MPI_Recv( evalk(:, ik), nwf, MPI_DOUBLE, irank, irank*nk+ik, MPI_COMM(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
        call MPI_Recv( Umnk(:, :, ik), nwf*eph_nwf_tot, MPI_DOUBLE_COMPLEX, irank, (mpi%procs+irank)*nk+ik, MPI_COMM(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
      end do
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
    use dfpt_variables, only: dfpt_kset, dfpt_Gkset, feveck0
    use dfpt_eigensystem, only: dfpt_eig_getevec
    use mod_wannier_variables, only: wf_fst, wf_lst, wf_nst, wf_nwf, wf_kset, wf_transform, wf_centers
    use mod_wannier_filehandling, only: wffile_readtransform
    use mod_wannier_helper, only: wfhelp_geteval, wfhelp_getevec
    use mod_kpointset, only: k_set, Gk_set
    use mod_eigensystem, only: nmatmax_ptr
    use mod_eigenvalue_occupancy, only: nstsv
    use mod_spin, only: nspinor
    use mod_APW_LO, only: nlotot
    use m_linalg, only: zlsp
    use xlapack, only: svd_divide_conquer
    use sorting, only: sort_index_1d
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

    integer, allocatable :: sort(:)
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
    ! (Resort energies. Might not be sorted in case of GW.)
    if (allocated(evalk)) deallocate( evalk )
    allocate( evalk(wf_fst:wf_lst, kset%nkpt) )
    allocate( sort(wf_fst:wf_lst) )
    do ik = 1, kset%nkpt
      call findkptinset( kset%vkl(:, ik), wf_kset, isym, ikw )
      call terminate_if_false( isym == 1, '(eph_el_read_wannier) &
        Requested k-point not found in Wannier k-point set.' )
      sort = sort_index_1d( wf_nst, eval(wf_fst:wf_lst, ikw) ) + wf_fst - 1
      evalk(:, ik) = eval(sort, ikw)
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
      allocate( evec_wan(nmatmax_ptr, nstsv, nspinor) )

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
  subroutine eph_el_gen_Hk_wannier( evalk, Umnk, Hk )
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

    CALL_ASSERT( is_square( Hk ),  '`Hk` must be square.' )
    CALL_ASSERT( nst == size( evalk ),  'Number of electronic states in `Umnk` and `evalk` do not match.' )
    CALL_ASSERT( nwf == size( Hk, dim=1 ),  'Number of Wannier functions in `Umnk` and `Hk` do not match.' )

    allocate( auxmat(nst, nwf) )
    do ist = 1, nst
      auxmat(ist, :) = evalk(ist) * Umnk(ist, :)
    end do
    call zgemm( 'c', 'n', nwf, nwf, nst, zone, Umnk, nst, auxmat, nst, zzero, Hk, nwf )
    deallocate( auxmat )
  end subroutine eph_el_gen_Hk_wannier

  !> Shifts a set of energies relative to the Fermi level and applies a scissor operator.
  pure subroutine eph_el_fermi_and_scissor( energies, ne, efermi, scissor )
    !> set of electron energies
    real(dp), intent(inout) :: energies(*)
    !> number of energies
    integer, intent(in) :: ne
    !> Fermi energy
    real(dp), intent(in) :: efermi
    !> scissor shift for occupied and unoccupied bands
    real(dp), intent(in) :: scissor(2)
  
    energies(:ne) = energies(:ne) - efermi
    if (abs(scissor(1)) < eph_el_degtol) then
      where (energies(:ne) > eph_el_degtol) 
        energies(:ne) = energies(:ne) + scissor(2)
      end where
    else if (abs(scissor(2)) < eph_el_degtol) then
      where (energies(:ne) < -eph_el_degtol) 
        energies(:ne) = energies(:ne) - scissor(1)
      end where
    else
      where (energies(:ne) < 0.0_dp) 
        energies(:ne) = energies(:ne) - scissor(1)
      end where
      where (energies(:ne) > 0.0_dp) 
        energies(:ne) = energies(:ne) + scissor(2)
      end where
    end if
  end subroutine eph_el_fermi_and_scissor

  !> Set default parameters for frequency grid used for electron energies.
  pure subroutine eph_el_set_default_frequency_grid( fgrid )
    use modinput, only: freq_grid_type
    !> frequency grid object
    type(freq_grid_type), intent(out) :: fgrid
  
    fgrid%type = 'density'          ! density based sampling
    fgrid%numpoints = 500           ! number of sampling points
    fgrid%range = [minval(eph_el_energy_k), maxval(eph_el_energy_k)]
    fgrid%padding = 0.1_dp          ! padding to add at both ends of range
    fgrid%lorentzwidth = 0.02_dp    ! width of Lorentzian density
  end subroutine eph_el_set_default_frequency_grid
  !-------------------------------------------------------------------------------- 
end module eph_electrons
