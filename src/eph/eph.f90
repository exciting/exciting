!> This is the main module for organizing and launching an
!> electron-phonon calculation using Wannier functions.
module eph
  use eph_variables

  use precision, only: dp
  use modmpi

  implicit none
  private

  public :: eph_launcher

contains

  !================================================================================ 
  ! MAIN TASK LAUNCHER
  !
  !> Set up and launch tasks of electron-phonon calculation.
  subroutine eph_launcher
    use eph_electrons
    use eph_phonons
    use eph_ephmat

    use dfpt_density_potential, only: dfpt_rhopot_init, dfpt_rhopot_free
    use dfpt_eigensystem, only: dfpt_eig_init, dfpt_eig_free
    use dfpt, only: dfpt_prepare, dfpt_finalize

    use modinput

    integer, parameter :: MAX_NUM_TASKS = 16

    integer :: num_tasks

    character(64) :: task_list(MAX_NUM_TASKS)

    !******************************************************************************** 
    ! CREATE LIST OF TASKS
    !
    num_tasks = 0

    ! compute EPH matrix elements on k- and q-grid in canonical coordinates
    if (associated(input%eph%ephmat)) then
      if (input%eph%ephmat%do /= 'skip') then
        num_tasks = num_tasks + 1 
        call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
          Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
        task_list(num_tasks) = 'gen_ephmat'
      end if
    end if

    ! interpolate electron energies on target set
    if (associated(input%eph%el_interpolation)) then
      num_tasks = num_tasks + 2
      call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
        Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
      task_list(num_tasks-1) = 'setup_interpolation_electrons'
      task_list(num_tasks) = 'interpolate_electrons'
    end if

    ! interpolate phonon energies on target set
    if (associated(input%eph%ph_interpolation)) then
      num_tasks = num_tasks + 2
      call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
        Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
      task_list(num_tasks-1) = 'setup_interpolation_phonons'
      task_list(num_tasks) = 'interpolate_phonons'
    end if

    ! interpolate EPH matrix to target set
    if (associated(input%eph%eph_interpolation)) then
      num_tasks = num_tasks + 2
      call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
        Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
      task_list(num_tasks-1) = 'setup_interpolation'
      task_list(num_tasks) = 'interpolate_ephmat'
    end if

    ! calculate electron self-energy Sigma and spectral function
    if (associated(input%eph%el_self_energy)) then
      if (input%eph%el_self_energy%do /= 'skip') then
        if (input%eph%el_self_energy%do == 'fromscratch') then
          num_tasks = num_tasks + 2
          call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
            Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
          task_list(num_tasks-1) = 'setup_interpolation'
          task_list(num_tasks) = 'gen_electron_selfenergy'
        end if
        num_tasks = num_tasks + 1 
        call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
          Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
        task_list(num_tasks) = 'write_electron_selfenergy_and_spectralfun'
      end if
    end if
    !................................................................................ 

    !******************************************************************************** 
    ! EXECUTE TASK LIST
    !
    ! initialize global EPH and DFPT variables
    call eph_var_init

    ! compute eigenvalues and eigenvectors on k-grid
    if (any( task_list == 'gen_ephmat' ) .or. &
        any( task_list == 'setup_interpolation' ) .or. &
        any( task_list == 'setup_interpolation_electrons' ) .or. &
        any( task_list == 'setup_interpolation_ephmat' )) then
      call dfpt_rhopot_init
      call dfpt_eig_init
      call dfpt_prepare
      call eph_el_set_energies
    end if

    ! generate EPH matrix elements on coarse grids
    if (any( task_list == 'gen_ephmat' )) &
      call eph_ephmat_gen_coarse( standard=input%eph%ephmat%standard )

    ! set up electron, phonon and EPH interpolation
    if (any( task_list == 'setup_interpolation' ) .or. &
        any( task_list == 'setup_interpolation_electrons' ) .or. &
        any( task_list == 'setup_interpolation_phonons' ) .or. &
        any( task_list == 'setup_interpolation_ephmat' )) &
      call eph_setup_interpolation( &
        electrons=any( task_list == 'setup_interpolation' ) .or. any( task_list == 'setup_interpolation_electrons' ), &
        phonons=any( task_list == 'setup_interpolation' ) .or. any( task_list == 'setup_interpolation_phonons' ), &
        ephmat=any( task_list == 'setup_interpolation' ) .or. any( task_list == 'setup_interpolation_ephmat' ) )

    ! interpolate electron energies
    if (any( task_list == 'interpolate_electrons' )) &
      call eph_interpolate_electrons

    ! interpolate phonon energies
    if (any( task_list == 'interpolate_phonons' )) &
      call eph_interpolate_phonons

    ! interpolate EPH matrix elements
    if (any( task_list == 'interpolate_ephmat' )) &
      call eph_interpolate_ephmat( include_polar=input%eph%eph_interpolation%include_polar )

    ! calculate electron self-energy Sigma
    if (any( task_list == 'gen_electron_selfenergy' )) &
      call eph_gen_electron_selfenergy

    ! write electron self-energy and spectral function to text file
    if (any( task_list == 'write_electron_selfenergy_and_spectralfun' )) &
      call eph_write_electron_selfenergy_specfun

    ! delete global EPH and DFPT variables
    if (any( task_list == 'gen_ephmat' ) .or. &
        any( task_list == 'setup_interpolation_electrons' )) then
      call dfpt_finalize( cleanup=.false. )
      call dfpt_eig_free
      call dfpt_rhopot_free
    end if
    call eph_el_free
    call eph_ph_free
    call eph_var_free
    !................................................................................ 
  end subroutine eph_launcher
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! SETUP WANNIER-FOURIER INTERPOLATION
  !
  !> Generate real-space matrices in localized atomic Wannier gauge for Wannier-Fourier-interpolation.
  subroutine eph_setup_interpolation( electrons, phonons, ephmat )
    use eph_electrons, only: eph_el_setup_interpolation
    use eph_phonons, only: eph_ph_setup_interpolation
    use eph_ephmat, only: eph_ephmat_setup_interpolation
    !> setup electron interpolation (default: `.true.`)
    logical, optional, intent(in) :: electrons
    !> setup phonons interpolation (default: `.true.`)
    logical, optional, intent(in) :: phonons
    !> setup EPH matrix interpolation (default: `.true.`)
    logical, optional, intent(in) :: ephmat
  
    logical :: el, ph, eph

    el = .true.
    if (present(electrons)) el = electrons
    ph = .true.
    if (present(phonons)) ph = phonons
    eph = .true.
    if (present(ephmat)) eph = ephmat

    ! electron interpolation
    if (el .or. eph) &
      call eph_el_setup_interpolation( write_localization=.true. )

    ! phonon interpolation
    if (ph .or. eph) &
      call eph_ph_setup_interpolation( write_localization=.true. )

    ! EPH matrix interpolation
    if (eph) &
      call eph_ephmat_setup_interpolation( write_localization=.true. )
  end subroutine eph_setup_interpolation
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! INTERPOLATE ELECTRON ENERGIES ON TARGET SET
  !
  !> Interpolate electron energies on target set and write result to file.
  subroutine eph_interpolate_electrons
    use eph_electrons, only: eph_el_interpolate
    use eph_inout, only: eph_io_write_energies

    use modinput

    real(dp), allocatable :: elengy(:,:)
    complex(dp), allocatable :: elevec(:,:,:)

    ! interpolate electron energies
    call eph_el_interpolate( eph_pset%vkl, elengy, elevec )
    ! write result to file
    if (eph_pset_path%num_points == 0) then
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, elengy, 'eph_el_energies', input%eph%el_interpolation%format, plist=eph_pset%vkl )
    else
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, elengy, 'eph_el_disp', input%eph%el_interpolation%format, path=eph_pset_path )
    end if
    if (allocated(elengy)) deallocate( elengy )
    if (allocated(elevec)) deallocate( elevec )
  end subroutine eph_interpolate_electrons
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! INTERPOLATE PHONON FREQUENCIES ON TARGET SET
  !
  !> Interpolate phonon frequencies on target set and write result to file.
  subroutine eph_interpolate_phonons
    use eph_phonons, only: eph_ph_interpolate
    use eph_inout, only: eph_io_write_energies

    use modinput

    real(kind=dp), allocatable :: phengy(:,:)
    complex(kind=dp), allocatable :: phevec(:,:,:)

    ! interpolate phonon frequencies
    call eph_ph_interpolate( eph_pset%vkl, phengy, phevec )
    ! write result to file
    if (eph_pset_path%num_points == 0) then
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, phengy, 'eph_ph_energies', input%eph%ph_interpolation%format, plist=eph_pset%vkl )
    else
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, phengy, 'eph_ph_disp', input%eph%ph_interpolation%format, path=eph_pset_path )
    end if
    if (allocated(phengy)) deallocate( phengy )
    if (allocated(phevec)) deallocate( phevec )
  end subroutine eph_interpolate_phonons
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! INTERPOLATE EPH MATRIX ON TARGET SET
  !
  !> Interpolate electron-phonon matrix elements on target set and write result to file.
  subroutine eph_interpolate_ephmat( include_polar )
    use eph_electrons, only: eph_el_interpolate
    use eph_phonons, only: eph_ph_interpolate
    use eph_ephmat, only: eph_ephmat_interpolate, eph_ephmat_average
    use eph_inout, only: eph_io_write_ephmat
  
    use modinput
    !> include (long-range) polar Fröhlich coupling in polar materials (default: `.true.`)
    logical, optional, intent(in) :: include_polar

    integer :: ik, iq, ikq, np, nk, nq, nst, nmode

    real(dp), allocatable, target :: vkl(:,:), vql(:,:), vkql(:,:), elengyk(:,:), elengykq(:,:), phengyq(:,:), gavg(:,:,:,:,:), gtmp(:,:,:)
    complex(dp), allocatable, target :: eleveck(:,:,:), eleveckq(:,:,:), phevecq(:,:,:), g(:,:,:,:,:)
    real(dp), pointer :: rptr1(:,:,:), rptr2(:,:,:)
    complex(dp), pointer :: zptr(:,:,:,:)

    ! set interpolation points
    np = eph_pset%nkpt
    if (input%eph%eph_interpolation%fix == 'k') then ! fixed k
      allocate( vkl(3, 1), vkql(3, np), vql(3, np) )
      vkl(:, 1) = input%eph%eph_interpolation%vplfix
      vql = eph_pset%vkl(:, :eph_pset%nkpt)
      vkql(1, :) = vql(1, :) + vkl(1, 1)
      vkql(2, :) = vql(2, :) + vkl(2, 1)
      vkql(3, :) = vql(3, :) + vkl(3, 1)
    else ! fixed q
      allocate( vkl(3, np), vkql(3, np), vql(3, 1) )
      vql(:, 1) = input%eph%eph_interpolation%vplfix
      vkl = eph_pset%vkl(:, :eph_pset%nkpt)
      vkql(1, :) = vkl(1, :) + vql(1, 1)
      vkql(2, :) = vkl(2, :) + vql(2, 1)
      vkql(3, :) = vkl(3, :) + vql(3, 1)
    end if

    ! interpolate electrons
    nk = size( vkl, dim=2 )
    call eph_el_interpolate( vkl, elengyk, eleveck )
    call eph_el_interpolate( vkql, elengykq, eleveckq )
    nst = size( elengyk, dim=1 )

    ! interpolate phonons
    nq = size( vql, dim=2 )
    call eph_ph_interpolate( vql, phengyq, phevecq )
    nmode = size( phengyq, dim=1 )

    ! interpolate EPH matrix
    zptr(1:size(eleveckq, dim=1), 1:size(eleveckq, dim=2), 1:nk, 1:nq) => eleveckq
    call eph_ephmat_interpolate( vkl, vql, eleveck, zptr, phengyq, phevecq, g, &
      include_polar=include_polar )

    ! average EPH matrix over degenerate states
    allocate( gavg(nst, nst, nmode, nk, nq) )
    !$omp parallel default( shared ) private( ik, iq, ikq, gtmp )
    !$omp do collapse(2)
    do iq = 1, nq
      do ik = 1, nk
        ikq = (iq - 1)*nk + ik
        call eph_ephmat_average( elengyk(:, ik), elengykq(:, ikq), input%eph%eph_interpolation%epsdegel, &
          phengyq(:, iq), input%eph%eph_interpolation%epsdegph, &
          g(:, :, :, ik, iq), gtmp )
        gavg(:, :, :, ik, iq) = gtmp
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! write results to file
    rptr1(1:3, 1:nk, 1:nq) => vkql ! map k+q vectors to 3D array
    rptr2(1:nst, 1:nk, 1:nq) => elengykq ! map k+q energies to 3D array
    if (mpiglobal%rank == 0) call eph_io_write_ephmat( eph_pset%bvec, vkl, rptr1, vql, elengyk, rptr2, phengyq, g, gavg, 'eph_ephmat', input%eph%eph_interpolation%format )

    deallocate( vkl, vql, vkql )
    deallocate( elengyk, elengykq, eleveck, eleveckq, phengyq, phevecq, g, gavg )
  end subroutine eph_interpolate_ephmat
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! ELECTRON SELF-ENERGY
  !
  !> Calculate electron-phonon contribution to electron self-energy \(\Sigma_{n{\bf k}}(\omega, T)\).
  !>
  !> This routine is MPI parallelized over electron bands \(n\) and phonon modes \(\nu\).
  subroutine eph_gen_electron_selfenergy
    use eph_variables, only: kset => eph_pset ! make target points locally accessible as `kset`
    use eph_electrons, only: eph_el_interpolate
    use eph_phonons, only: eph_ph_interpolate
    use eph_ephmat, only: eph_ephmat_interpolate
    use eph_electron_selfenergy

    use block_data_file, only: block_data_file_type
    use grid_utils, only: linspace
    use mod_kpointset, only: k_set, delete_k_vectors
    use mod_opt_tetra, only: t_set, opt_tetra_destroy
    use mod_mpi_allreduce, only: xmpi_allreduce
    use modinput

    real(dp), parameter :: vgamma(3, 1) = 0.0_dp  ! Gamma point

    integer :: fmode, lmode, fwf, lwf, ntemp, nfreq, i, ik, itemp
    logical :: autofreq
    type(freq_grid_type) :: fgrid
    type(k_set) :: qset
    type(t_set) :: tset
    type(block_data_file_type) :: selfen_file
    character(:), allocatable :: fname

    integer, allocatable :: mode_indices(:,:), &  ! phonon mode indices in distributed calculation
                            band_indices(:,:)     ! band indices in distributed calculation
    real(dp), allocatable :: temps(:), &          ! temperatures
                             freqs(:), &          ! frequencies
                             vkql(:,:), &         ! shifted BZ integration set k+q
                             el_energy_k(:,:), &  ! electron energies on target set k
                             el_energy_kq(:,:), & ! electron energies on shifted BZ integration set k+q
                             ph_energy_q(:,:)     ! phonon frequencies on BZ integration set q
    complex(dp), allocatable, target :: el_evec_k(:,:,:), &   ! electron eigenvectors U(k) on target set k
                                        el_evec_kq(:,:,:), &  ! electron eigenvectors U(k+q) on shifted BZ integration set k+q
                                        ph_evec_q(:,:,:), &   ! phonon eigenvectors e(q) on BZ integration set q
                                        ephmat(:,:,:,:,:), &  ! EPH matrix elements g(k,q) on BZ integration set
                                        ephmat0(:,:,:,:,:), & ! EPH matrix elements g(k,0) in atomic gauge on target set k
                                        selfen_fm(:,:,:), &   ! Fan-Migdal self-energy
                                        selfen_dw(:,:)        ! Debye-Waller electron self-energy  
    complex(dp), pointer :: ptr(:,:,:,:)

    ! distribute phonon modes and electron bands among processes
    call patchwork_distribution( eph_nmode, eph_nwf, mpiglobal%procs, mode_indices, band_indices )
    fmode = eph_fmode + mode_indices(1, mpiglobal%rank+1) - 1
    lmode = eph_fmode + mode_indices(2, mpiglobal%rank+1) - 1
    fwf = eph_fwf + band_indices(1, mpiglobal%rank+1) - 1
    lwf = eph_fwf + band_indices(2, mpiglobal%rank+1) - 1
    ! set frequency grid
    call eph_else_set_default_frequency_grid( fgrid )
    if (associated(input%eph%el_self_energy%freq_grid)) &
      fgrid = input%eph%el_self_energy%freq_grid
    nfreq = fgrid%numpoints
    autofreq = all( fgrid%range == 0.0_dp )
    ! set temperatures
    temps = linspace( max( eph_temperature_zero, input%eph%el_self_energy%tempset(1) ), &
                      max( eph_temperature_zero, input%eph%el_self_energy%tempset(2) ), &
                      input%eph%el_self_energy%tempset(3) )
    ntemp = size(temps)
    ! initialize BZ integration points q
    call eph_var_init_bz_int( input%eph%el_self_energy%ngridbz, input%eph%el_self_energy%vbzoff, qset, tset )
    ! interpolate electrons on target points k
    call eph_el_interpolate( kset%vkl(:, 1:kset%nkpt), el_energy_k, el_evec_k )
    ! interpolate phonons on BZ integration points q
    call eph_ph_interpolate( qset%vkl(:, 1:qset%nkpt), ph_energy_q, ph_evec_q )
    ! set up and open binary file for self-energy
    fname = eph_else_get_binary_file_name( qset%ngridk, [eph_fwf, eph_lwf], nfreq, temps )
    selfen_file = block_data_file_type( fname, [nfreq+1, eph_nwf+1, ntemp], cmplx( 0, 0, dp ) )
    call selfen_file%open( mpiglobal, delete_existing=.true. )

    allocate( selfen_fm(nfreq+1, eph_nwf+1, ntemp), selfen_dw(eph_nwf, ntemp) )

    ! loop over target points k
    do ik = 1, kset%nkpt
      ! generate k+q vectors
      vkql = reshape( [(qset%vkl(:, i)+kset%vkl(:, ik), i=1, qset%nkpt)], [3, qset%nkpt] )
      ! interpolate electrons on k+q
      call eph_el_interpolate( vkql, el_energy_kq, el_evec_kq )
      ! generate frequency grid
      if (autofreq) fgrid%range = [minval( el_energy_kq ), maxval( el_energy_kq )]
      freqs = eph_var_gen_frequency_grid( fgrid, el_energy_k(:, ik) )

      !******************************************************************************** 
      ! Fan-Migdal self-energy
      !******************************************************************************** 
      selfen_fm = cmplx( 0, 0, dp )
      ! interpolate EPH matrix on (k,q)
      ptr(1:size(el_evec_kq, dim=1), 1:size(el_evec_kq, dim=2), 1:1, 1:size(el_evec_kq, dim=3)) => el_evec_kq
      call eph_ephmat_interpolate( kset%vkl(:, ik:ik), qset%vkl(:, 1:qset%nkpt), &
        el_evec_k(:, :, ik:ik), ptr, ph_energy_q, ph_evec_q, &
        ephmat, &
        frange=[fwf, lwf]-eph_fwf+1, mrange=[fmode, lmode]-eph_fmode+1 )
      ! compute self-energy
      ptr(1:size(ephmat, dim=1), 1:size(ephmat, dim=2), 1:size(ephmat, dim=3), 1:size(ephmat, dim=5)) => ephmat
      select case (input%eph%el_self_energy%integration)
        ! smearing
        case ('smearing')
          call eph_else_gen_fan_migdal_smearing( freqs, temps, el_energy_kq(fwf:lwf, :), ph_energy_q(fmode:lmode, :), &
            ptr, qset, input%eph%el_self_energy%swidth, selfen_fm  )
        ! Kramers-Kronig
        case default
          call eph_else_gen_fan_migdal_aux( freqs, temps, el_energy_kq(fwf:lwf, :), ph_energy_q(fmode:lmode, :), &
            ptr, tset, selfen_fm  )
      end select
      call xmpi_allreduce( selfen_fm, mpiglobal )

      !******************************************************************************** 
      ! Debye-Waller self-energy
      !******************************************************************************** 
      selfen_dw = cmplx( 0, 0, dp )
      ! interpolate EPH matrix on (k,0) in atomic gauge
      ptr(1:size(el_evec_k, dim=1), 1:size(el_evec_k, dim=2), 1:1, 1:1) => el_evec_k(:, :, ik)
      call eph_ephmat_interpolate( kset %vkl(:, ik:ik), vgamma, &
        el_evec_k(:, :, ik:ik), ptr, ph_energy_q(:, 1:1), ph_evec_q(:, :, 1:1), &
        ephmat0, &
        phonon_gauge='a', frange=[fwf, lwf]-eph_fwf+1 )
      ! compute self-energy
      call eph_else_gen_debye_waller_ahc( temps, el_energy_k(:, ik), el_energy_k(fwf:lwf, ik), ph_energy_q(fmode:lmode, :), ph_evec_q(:, fmode:lmode, :), &
        ephmat0(:, :, :, 1, 1), qset, input%eph%el_self_energy%swidth, selfen_dw )
      call xmpi_allreduce( selfen_dw, mpiglobal )

      ! store additional data (frequency grid, electron energies and DW self-energy)
      selfen_fm(:nfreq, eph_nwf+1, 1) = cmplx( freqs, 0, dp )
      do itemp = 1, ntemp
        selfen_fm(nfreq+1, :eph_nwf, itemp) = cmplx( el_energy_k(:, ik), selfen_dw(:, itemp)%re, dp )
      end do
      ! write electron self-energy to file
      call selfen_file%write( ik, selfen_fm )
    end do

    ! close binary file
    call selfen_file%close( mpiglobal )

    ! clean up
    call delete_k_vectors( qset )
    call opt_tetra_destroy( tset )
  end subroutine eph_gen_electron_selfenergy
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! ELECTRON SPECTRAL FUNCTION
  !
  !> Read electron self-energy from binary file and evaluate it and the corresponding spectral function
  !> on a different output frequency grid and write the result to a text file.
  subroutine eph_write_electron_selfenergy_specfun
    use eph_variables, only: kset => eph_pset ! make target points locally accessible as `kset`
    use eph_electron_selfenergy

    use grid_utils, only: linspace
    use block_data_file, only: block_data_file_type
    use os_utils, only: make_directory
    use convolution, only: smoothen
    use modinput

    integer :: nfreq_in, nfreq_out, ntemp, un, ik, ik1, ik2, itemp, ist, ifreq, stat
    logical :: autofreq
    type(freq_grid_type) :: fgrid_in, fgrid_out
    type(block_data_file_type) :: selfen_file
    character(:), allocatable :: fname
    character(128) :: fn, dn

    real(dp), allocatable :: temps(:), el_energy_k(:), selfen_dw(:), &
                             freqs_in(:), freqs_out(:, :), &
                             sfun_in(:, :), sfun_out(:, :), &
                             fun1(:), fun2(:)
    complex(dp), allocatable :: selfen_in(:,:,:), selfen_out(:,:)

    ! set input frequency grid
    call eph_else_set_default_frequency_grid( fgrid_in )
    if (associated(input%eph%el_self_energy%freq_grid)) &
      fgrid_in = input%eph%el_self_energy%freq_grid
    nfreq_in = fgrid_in%numpoints
    ! set output frequency grid
    call eph_else_set_default_frequency_grid( fgrid_out )
    if (associated(input%eph%el_self_energy%output_settings)) then 
      if (associated(input%eph%el_self_energy%output_settings%freq_grid)) &
        fgrid_out = input%eph%el_self_energy%output_settings%freq_grid
        fgrid_out%padding = 0.0_dp
    end if
    nfreq_out = fgrid_out%numpoints
    autofreq = all( fgrid_out%range == 0.0_dp )
    ! set temperatures
    temps = linspace( max( eph_temperature_zero, input%eph%el_self_energy%tempset(1) ), &
                      max( eph_temperature_zero, input%eph%el_self_energy%tempset(2) ), &
                      input%eph%el_self_energy%tempset(3) )
    ntemp = size(temps)
    ! set output directory
    dn = '.'
    if (associated(input%eph%el_self_energy%output_settings)) &
      dn = trim( adjustl( input%eph%el_self_energy%output_settings%directory ) )
    if (trim( adjustl( dn ) ) /= '.') then
      un = make_directory( dn, comm=mpiglobal )
      call terminate_if_false( mpiglobal, un==0, 'Failed to create output directory "'//trim( adjustl( dn ) )//'".' )
    end if
    ! open binary file with electron self-energy
    fname = eph_else_get_binary_file_name( input%eph%el_self_energy%ngridbz, [eph_fwf, eph_lwf], nfreq_in, temps )
    selfen_file = block_data_file_type( fname, [nfreq_in+1, eph_nwf+1, ntemp], cmplx( 0, 0, dp ) )
    call selfen_file%open( mpiglobal )
    ! set loop limits
    ik1 = firstofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
    ik2 = lastofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )

    allocate( selfen_in(nfreq_in+1, eph_nwf+1, ntemp), selfen_dw(eph_nwf) )
    allocate( freqs_out(nfreq_out, eph_nwf), selfen_out(nfreq_out, eph_nwf) )
    allocate( sfun_in(nfreq_in, eph_nwf), sfun_out(nfreq_out, eph_nwf) )

    do ik = ik1, ik2
      ! read input data from file
      call selfen_file%read( ik, selfen_in )
      freqs_in = selfen_in(:nfreq_in, eph_nwf+1, 1)%re
      el_energy_k = selfen_in(nfreq_in+1, :eph_nwf, 1)%re
      if (autofreq) fgrid_out%range = [freqs_in(1), freqs_in(nfreq_in)]

      do itemp = 1, ntemp
        selfen_dw = selfen_in(nfreq_in+1, :eph_nwf, itemp)%im

        do ist = 1, eph_nwf
          ! apply Lorentzian broadening for integration without smearing
          if (input%eph%el_self_energy%integration /= 'smearing') then
            fun1 = selfen_in(:nfreq_in, ist, itemp)%re
            fun2 = selfen_in(:nfreq_in, ist, itemp)%im
            call smoothen( freqs_in, fun1, input%eph%el_self_energy%swidth, kernel='lorentzian', mode='linear' )
            call smoothen( freqs_in, fun2, input%eph%el_self_energy%swidth, kernel='lorentzian', mode='linear' )
            selfen_in(:nfreq_in, ist, itemp) = cmplx( fun1, fun2, dp )
          end if
          ! obtain true self-energy from auxiliary one
          if (input%eph%el_self_energy%integration == 'kramers-kronig') then
            call eph_else_gen_fan_migdal_from_aux( freqs_in, selfen_in(:nfreq_in, ist, itemp) )
          end if
          ! add FM and DW self-energy
          selfen_in(:nfreq_in, ist, itemp) = selfen_in(:nfreq_in, ist, itemp) + selfen_dw(ist) 
          ! resample self-energy
          fgrid_out%range = fgrid_out%range - el_energy_k(ist)
          call eph_else_resample( nfreq_in, freqs_in-el_energy_k(ist), selfen_in(:nfreq_in, ist, itemp), fgrid_out, freqs_out(:, ist), selfen_out(:, ist), &
            interpolation_method=input%eph%el_self_energy%output_settings%interpolation_method )
          freqs_out(:, ist) = freqs_out(:, ist) + el_energy_k(ist) 
          fgrid_out%range = fgrid_out%range + el_energy_k(ist)
          ! compute spectral function
          call eph_else_gen_specfun( freqs_in-el_energy_k(ist), selfen_in(:nfreq_in, ist, itemp), sfun_in(:, ist) )
          call eph_else_gen_specfun( freqs_out(:, ist)-el_energy_k(ist), selfen_out(:, ist), sfun_out(:, ist) )
        end do

        ! write self-energy and spectral function to file
        write( fn, '("EPH_ELSE+SFUN_P",i3.3,"_T",i4.4,".OUT")' ) ik, nint( temps(itemp) )
        open( newunit=un, file=trim(dn)//'/'//trim(fn), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_write_electron_selfenergy_specfun) &
          Failed to open file `'//trim(dn)//'/'//trim(fn)//'`.' )
        write( un, '("#",10000g26.16e3)' ) el_energy_k
        write( un, '("#",10000g26.16e3)' ) selfen_dw
        do ifreq = 1, nfreq_in
          do ist = 1, eph_nwf
            write( un, '(4g26.16e3)', advance='no' ) freqs_in(ifreq), selfen_in(ifreq, ist, itemp), sfun_in(ifreq, ist)
          end do
          write( un, * )
        end do
        close( un )

        write( fn, '("EPH_ELSE+SFUN_RSMPLD_P",i3.3,"_T",i4.4,".OUT")' ) ik, nint( temps(itemp) )
        open( newunit=un, file=trim(dn)//'/'//trim(fn), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_write_electron_selfenergy_specfun) &
          Failed to open file `'//trim(dn)//'/'//trim(fn)//'`.' )
        write( un, '("#",10000g26.16e3)' ) el_energy_k
        write( un, '("#",10000g26.16e3)' ) selfen_dw
        do ifreq = 1, nfreq_out
          do ist = 1, eph_nwf
            write( un, '(4g26.16e3)', advance='no' ) freqs_out(ifreq, ist), selfen_out(ifreq, ist), sfun_out(ifreq, ist)
          end do
          write( un, * )
        end do
        close( un )

      end do
    end do

    deallocate( selfen_in, selfen_dw, selfen_out, freqs_out, sfun_in, sfun_out )

    ! close binary file
    call selfen_file%close( mpiglobal )
  end subroutine eph_write_electron_selfenergy_specfun
  !-------------------------------------------------------------------------------- 
end module eph
