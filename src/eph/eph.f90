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
        if (input%eph%ephmat%hilo) then ! for Wannier band map and phonon eigenvectors
          num_tasks = num_tasks + 2 
          call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
            Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
          task_list(num_tasks-1) = 'setup_interpolation_electrons'
          task_list(num_tasks) = 'setup_interpolation_phonons'
        end if
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

    ! calculate electron self-energy Sigma, QP energies and spectral function
    if (associated(input%eph%el_self_energy)) then
      if (input%eph%el_self_energy%do /= 'skip') then
        if (input%eph%el_self_energy%do == 'fromscratch') then
          num_tasks = num_tasks + 2
          call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
            Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
          task_list(num_tasks-1) = 'setup_interpolation'
          task_list(num_tasks) = 'gen_electron_selfenergy'
        end if
        num_tasks = num_tasks + 2 
        call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
          Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
        task_list(num_tasks-1) = 'write_electron_selfenergy_and_spectralfun'
        task_list(num_tasks) = 'solve_electron_qp_equation'
      end if
    end if

    ! calculate EPH coupling strength for electrons
    if (associated(input%eph%coupling_strength)) then
      if (input%eph%coupling_strength%particle == 'electron' .and. input%eph%coupling_strength%do == 'fromscratch') then
        num_tasks = num_tasks + 2
        call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
          Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
        task_list(num_tasks-1) = 'setup_interpolation'
        task_list(num_tasks) = 'gen_electron_coupling_strength'
      end if
    end if

    ! calculate EPH coupling strength for phonons
    if (associated(input%eph%coupling_strength)) then
      if (input%eph%coupling_strength%particle == 'phonon' .and. input%eph%coupling_strength%do == 'fromscratch') then
        num_tasks = num_tasks + 2
        call terminate_if_false( num_tasks <= MAX_NUM_TASKS, '(eph_launcher) &
          Maximum number of tasks exceeded. Increase `MAX_NUM_TASKS`.' )
        task_list(num_tasks-1) = 'setup_interpolation'
        task_list(num_tasks) = 'gen_phonon_coupling_strength'
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
      call eph_el_set_wannier_eigensystem
    end if

    ! set up electron interpolation
    if (any( task_list == 'setup_interpolation' ) .or. &
        any( task_list == 'setup_interpolation_electrons' )) &
      call eph_setup_interpolation( electrons=.true. )

    ! set up phonon interpolation
    if (any( task_list == 'setup_interpolation' ) .or. &
        any( task_list == 'setup_interpolation_phonons' )) &
      call eph_setup_interpolation( phonons=.true. )

    ! generate EPH matrix elements on coarse grids
    if (any( task_list == 'gen_ephmat' )) &
      call eph_ephmat_gen_coarse( standard=input%eph%ephmat%standard, hilo=input%eph%ephmat%hilo )

    ! set up EPH interpolation
    if (any( task_list == 'setup_interpolation' ) .or. &
        any( task_list == 'setup_interpolation_ephmat' )) &
      call eph_setup_interpolation( ephmat=.true. )

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

    ! solve electron quasi-particle equation
    if (any( task_list == 'solve_electron_qp_equation' )) &
      call eph_solve_electron_quasi_particle_equation

    ! calculate electron EPH coupling strength
    if (any( task_list == 'gen_electron_coupling_strength' )) &
      call eph_gen_electron_coupling_strength

    ! calculate phonon EPH coupling strength
    if (any( task_list == 'gen_phonon_coupling_strength' )) &
      call eph_gen_phonon_coupling_strength

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
    !> setup electron interpolation (default: `.false.`)
    logical, optional, intent(in) :: electrons
    !> setup phonons interpolation (default: `.false.`)
    logical, optional, intent(in) :: phonons
    !> setup EPH matrix interpolation (default: `.false.`)
    logical, optional, intent(in) :: ephmat
  
    logical :: el, ph, eph

    el = .false.
    if (present(electrons)) el = electrons
    ph = .false.
    if (present(phonons)) ph = phonons
    eph = .false.
    if (present(ephmat)) eph = ephmat

    ! electron interpolation
    if (el) &
      call eph_el_setup_interpolation( write_localization=.true. )

    ! phonon interpolation
    if (ph) &
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
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, eph_fst_span, elengy, 'eph_el_energies', input%eph%el_interpolation%format, plist=eph_pset%vkl(:, 1:eph_pset%nkpt) )
    else
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, eph_fst_span, elengy, 'eph_el_disp', input%eph%el_interpolation%format, path=eph_pset_path )
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
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, eph_fmode, phengy, 'eph_ph_energies', input%eph%ph_interpolation%format, plist=eph_pset%vkl(:, 1:eph_pset%nkpt) )
    else
      if (mpiglobal%rank == 0) call eph_io_write_energies( eph_pset%bvec, eph_fmode, phengy, 'eph_ph_disp', input%eph%ph_interpolation%format, path=eph_pset_path )
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
        call eph_ephmat_average( elengyk(:, ik), elengykq(:, ikq), eph_el_degtol, &
          phengyq(:, iq), eph_ph_degtol, &
          g(:, :, :, ik, iq), gtmp )
        gavg(:, :, :, ik, iq) = gtmp
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! write results to file
    rptr1(1:3, 1:nk, 1:nq) => vkql ! map k+q vectors to 3D array
    rptr2(1:nst, 1:nk, 1:nq) => elengykq ! map k+q energies to 3D array
    if (mpiglobal%rank == 0) call eph_io_write_ephmat( eph_pset%bvec, vkl, rptr1, vql, eph_fst_span, elengyk, rptr2, eph_fmode, phengyq, g, gavg, 'eph_ephmat', input%eph%eph_interpolation%format )

    deallocate( vkl, vql, vkql )
    deallocate( elengyk, elengykq, eleveck, eleveckq, phengyq, phevecq, g, gavg )
  end subroutine eph_interpolate_ephmat
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! ELECTRON SELF-ENERGY
  !
  !> Calculate electron-phonon contribution to electron self-energy \(\Sigma_{n{\bf k}}(\omega, T)\).
  !>
  !> This routine is MPI parallelized over electron bands \(m\) and phonon modes \(\nu\).
  subroutine eph_gen_electron_selfenergy
    use eph_variables, only: kset => eph_pset ! make target points locally accessible as `kset`
    use eph_electrons, only: eph_el_interpolate
    use eph_phonons, only: eph_ph_interpolate
    use eph_ephmat, only: eph_ephmat_interpolate
    use eph_inout, only: eph_io_write_el_self_energy
    use eph_electron_selfenergy

    use block_data_file, only: block_data_file_type
    use grid_utils, only: linspace
    use mod_kpointset, only: k_set, generate_k_vectors, delete_k_vectors
    use mod_opt_tetra, only: t_set, opt_tetra_destroy
    use mod_mpi_allreduce, only: xmpi_allreduce
    use modinput

    real(dp), parameter :: vgamma(3, 1) = 0.0_dp  ! Gamma point

    integer :: fmode, lmode, fwf, lwf, ntemp, nfreq, i, ik, itemp
    real(dp) :: swidth
    type(freq_grid_type) :: fgrid
    type(k_set) :: qset
    type(t_set) :: tset
    type(block_data_file_type) :: selfen_file

    integer, allocatable :: mode_indices(:,:), &  ! phonon mode indices in distributed calculation
                            band_indices(:,:)     ! band indices in distributed calculation
    real(dp), allocatable :: temps(:), &          ! temperatures
                             freqs(:), &          ! frequencies
                             dfreqs(:), &         ! base frequency sampling density
                             vkql(:,:), &         ! shifted BZ integration set k+q
                             el_energy_k0(:,:), & ! electron energies on coarse original k-grid
                             el_energy_k(:,:), &  ! electron energies on target set k
                             el_energy_kq(:,:), & ! electron energies on shifted BZ integration set k+q
                             ph_energy_q(:,:)     ! phonon frequencies on BZ integration set q
    complex(dp), allocatable, target :: el_evec_k(:,:,:), &         ! electron eigenvectors U(k) on target set k
                                        el_evec_kq(:,:,:), &        ! electron eigenvectors U(k+q) on shifted BZ integration set k+q
                                        ph_evec_q(:,:,:), &         ! phonon eigenvectors e(q) on BZ integration set q
                                        ephmat(:,:,:,:,:), &        ! EPH matrix elements g(k,q) on BZ integration set
                                        ephmat0(:,:,:,:,:), &       ! EPH matrix elements g(k,0) in atomic gauge on target set k
                                        selfen_fm(:,:,:), &         ! Fan-Migdal self-energy
                                        selfen_dw(:,:), &           ! Debye-Waller electron self-energy  
                                        selfen_fm_hilo(:,:,:), &    ! high and low energy Fan-Migdal electron self-energy
                                        selfen_dw_hilo(:,:,:), &    ! high and low energy Debye-Waller electron self-energy
                                        selfen_fm_hilo_R(:,:,:,:), &! high and low energy Fan-Migdal electron self-energy in real-space Wannier gauge
                                        selfen_dw_hilo_R(:,:,:,:)   ! high and low energy Debye-Waller electron self-energy in real-space Wannier gauge
    complex(dp), pointer :: ptr(:,:,:,:)

    ! set target k-points
    if (associated(input%eph%target)) then
      call eph_var_read_target
    else
      call generate_k_vectors( eph_pset, eph_kset_el%bvec, eph_kset_el%ngridk, eph_kset_el%vkloff, .true., uselibzint=.false. )
    end if
    ! distribute phonon modes and electron bands among processes
    call patchwork_distribution( eph_nmode, eph_nwf, mpiglobal%procs, mode_indices, band_indices )
    fmode = mode_indices(1, mpiglobal%rank+1) + eph_fmode - 1
    lmode = mode_indices(2, mpiglobal%rank+1) + eph_fmode - 1
    fwf = band_indices(1, mpiglobal%rank+1) + eph_fwf - 1
    lwf = band_indices(2, mpiglobal%rank+1) + eph_fwf - 1
    ! set temperatures
    temps = linspace( max( eph_temperature_zero, input%eph%el_self_energy%tempset(1) ), &
                      max( eph_temperature_zero, input%eph%el_self_energy%tempset(2) ), &
                      input%eph%el_self_energy%tempset(3) )
    ntemp = size(temps)
    ! initialize BZ integration points q
    call eph_var_init_bz_int( input%eph%el_self_energy%ngridbz, input%eph%el_self_energy%vbzoff, .false., qset, tset )
    ! interpolate electrons on coarse original k-grid
    call eph_el_interpolate( eph_kset_el%vkl(:, 1:eph_kset_el%nkpt), el_energy_k0, el_evec_k )
    ! interpolate electrons on target points k
    call eph_el_interpolate( kset%vkl(:, 1:kset%nkpt), el_energy_k, el_evec_k )
    ! interpolate phonons on BZ integration points q
    call eph_ph_interpolate( qset%vkl(:, 1:qset%nkpt), ph_energy_q, ph_evec_q, mrange=[fmode, lmode] )
    ! interpolate high and low energy self-energy on target points k
    call eph_else_hilo_setup_interpolation( temps, selfen_fm_hilo_R, selfen_dw_hilo_R )
    call eph_else_hilo_interpolate( kset%vkl(:, 1:kset%nkpt), el_energy_k, el_evec_k, selfen_fm_hilo_R, selfen_fm_hilo )
    call eph_else_hilo_interpolate( kset%vkl(:, 1:kset%nkpt), el_energy_k, el_evec_k, selfen_dw_hilo_R, selfen_dw_hilo )
    ! set frequency grid
    ! Note: The grid should have the following properties:
    !       1. no sampling where the DOS is zero
    !       2. densely sampled around the e_nk of the current k-point
    !       3. also sufficiently dense elsewhere to be able to apply the Hilbert transform
    !          from imaginary to real part (Hilbert transform is a non-local operation and needs the full frequency space)
    !       This is achieved by generating a density based sampling where one part of the density is generated from
    !       the electron energies on the original coarse k-grid in the full BZ (for properties 1 and 3) and another
    !       part of the density is based on the electron energies of the current target k-point (property 2).
    call eph_else_set_default_frequency_grid( fgrid )
    fgrid%numpoints = input%eph%el_self_energy%nfreqperband * eph_nwf
    if (associated(input%eph%el_self_energy%freq_grid)) &
      fgrid = input%eph%el_self_energy%freq_grid
    nfreq = fgrid%numpoints
    if (all( fgrid%range == 0.0_dp )) fgrid%range = [minval( el_energy_k0 ), maxval( el_energy_k0 )]
    freqs = eph_var_gen_frequency_grid( fgrid, reshape( el_energy_k0, [size( el_energy_k0 )] ) )
    allocate( dfreqs(nfreq) )
    dfreqs(2:nfreq-1) = 2.0_dp / (freqs(3:) - freqs(:nfreq-2))
    dfreqs(1) = 1.0_dp / (freqs(2) - freqs(1))
    dfreqs(nfreq) = 1.0_dp / (freqs(nfreq) - freqs(nfreq-1))
    ! set up and open binary file for self-energy
    selfen_file = block_data_file_type( eph_else_filename, [nfreq+2, eph_nwf+1, ntemp], cmplx( 0, 0, dp ) )
    call selfen_file%open( mpiglobal, delete_existing=.true. )

    allocate( selfen_fm(nfreq, eph_fwf:eph_lwf, ntemp), selfen_dw(eph_fwf:eph_lwf, ntemp) )

    ! loop over target points k
    do ik = 1, kset%nkpt
      ! generate frequency grid
      freqs = eph_var_gen_frequency_grid( fgrid, energies=el_energy_k(:, ik), density=dfreqs )
      ! generate k+q vectors
      vkql = reshape( [(qset%vkl(:, i)+kset%vkl(:, ik), i=1, qset%nkpt)], [3, qset%nkpt] )
      ! interpolate electrons on k+q
      call eph_el_interpolate( vkql, el_energy_kq, el_evec_kq, irange=[fwf, lwf] )

      !******************************************************************************** 
      ! Fan-Migdal self-energy
      !******************************************************************************** 
      selfen_fm = cmplx( 0, 0, dp )
      ! interpolate EPH matrix on (k,q)
      ptr(1:size(el_evec_kq, dim=1), 1:size(el_evec_kq, dim=2), 1:1, 1:size(el_evec_kq, dim=3)) => el_evec_kq
      call eph_ephmat_interpolate( kset%vkl(:, ik:ik), qset%vkl(:, 1:qset%nkpt), &
        el_evec_k(:, :, ik:ik), ptr, ph_energy_q, ph_evec_q, &
        ephmat, &
        frange=[fwf, lwf], mrange=[fmode, lmode] )
      ! compute self-energy
      ptr(1:size(ephmat, dim=1), 1:size(ephmat, dim=2), 1:size(ephmat, dim=3), 1:size(ephmat, dim=5)) => ephmat
      select case (input%eph%el_self_energy%integration)
        ! smearing
        case ('smearing')
          swidth = input%eph%el_self_energy%swidth
          call eph_else_gen_fan_migdal_smearing( freqs, temps, el_energy_k(:, ik), el_energy_kq, ph_energy_q, &
            ptr, qset, swidth, selfen_fm )
        ! Kramers-Kronig
        case ('kramers-kronig')
          swidth = 0.0_dp
          call eph_else_gen_fan_migdal_aux( freqs, temps, el_energy_k(:, ik), el_energy_kq, ph_energy_q, &
            ptr, tset, selfen_fm )
        case default
          call terminate_if_false( .false., '(eph_gen_electron_selfenergy) &
            Invalid integration method.' )
      end select
      call xmpi_allreduce( selfen_fm, mpiglobal )

      !******************************************************************************** 
      ! Debye-Waller self-energy
      !******************************************************************************** 
      selfen_dw = cmplx( 0, 0, dp )
      ! interpolate EPH matrix on (k,0) in atomic gauge
	    ptr(lbound(el_evec_k, dim=1):ubound(el_evec_k, dim=1), 1:size(el_evec_k, dim=2), 1:1, 1:1) => el_evec_k(:, :, ik)
      call eph_ephmat_interpolate( kset%vkl(:, ik:ik), vgamma, &
        el_evec_k(:, :, ik:ik), ptr(fwf:lwf, :, :, :), ph_energy_q(:, 1:1), ph_evec_q(:, :, 1:1), &
        ephmat0, &
        phonon_gauge='a', frange=[fwf, lwf] )
      ! compute self-energy
      call eph_else_gen_debye_waller_ahc( temps, el_energy_k(:, ik), el_energy_k(fwf:lwf, ik), ph_energy_q, ph_evec_q, &
        ephmat0(:, :, :, 1, 1), qset, input%eph%el_self_energy%swidth, selfen_dw )
      call xmpi_allreduce( selfen_dw, mpiglobal )

      ! write electron self-energy to file
      call eph_io_write_el_self_energy( selfen_file, ik, eph_fst_span, freqs, temps, el_energy_k(:, ik), selfen_fm, selfen_dw%re, selfen_fm_hilo(:, :, ik)%re, selfen_dw_hilo(:, :, ik)%re, &
        input%eph%el_self_energy%integration, swidth )
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
    use eph_inout, only: eph_io_write_else_sfun
    use eph_electron_selfenergy

    use grid_utils, only: linspace
    use block_data_file, only: block_data_file_type
    use os_utils, only: make_directory_pure
    use mod_kpointset, only: generate_k_vectors
    use modinput

    integer :: nfreq_in, nfreq_out, ntemp, fst, lst, ik, ik1, ik2, itemp, ist, stat 
    logical :: autofreq
    type(freq_grid_type) :: fgrid_out
    type(block_data_file_type) :: selfen_file
    character(:), allocatable :: string, interpolation_method, format
    character(128) :: fname, dirname

    real(dp), allocatable :: temps(:), el_energy_k(:), selfen_dw(:,:), &
                             freqs_in(:), freqs_out(:,:), &
                             sfun_in(:,:), sfun_out(:,:)
    complex(dp), allocatable :: selfen_in(:,:,:), selfen_out(:,:)

    ! set target k-points
    if (associated(input%eph%target)) then
      call eph_var_read_target
    else
      call generate_k_vectors( eph_pset, eph_kset_el%bvec, eph_kset_el%ngridk, eph_kset_el%vkloff, .true., uselibzint=.false. )
    end if
    ! open binary file with electron self-energy
    selfen_file = block_data_file_type( eph_else_filename, [-1], cmplx( 0, 0, dp ) )
    call selfen_file%open( mpiglobal )
    ! read input frequency grid and temperatures from file
    call eph_else_gen_from_file( selfen_file, 1, freqs_in, temps, el_energy_k, selfen_in, selfen_dw )
    nfreq_in = size( freqs_in )
    ntemp = size( temps )
    fst = lbound( el_energy_k, dim=1 )
    lst = ubound( el_energy_k, dim=1 )
    ! set output frequency grid
    call eph_else_set_default_frequency_grid( fgrid_out )
    if (associated(input%eph%el_self_energy%output_settings)) then 
      if (associated(input%eph%el_self_energy%output_settings%freq_grid)) &
        fgrid_out = input%eph%el_self_energy%output_settings%freq_grid
        fgrid_out%padding = 0.0_dp
    end if
    nfreq_out = fgrid_out%numpoints
    autofreq = all( fgrid_out%range == 0.0_dp )
    ! set output settings
    string = eph_else_setting_string( input%eph%el_self_energy%ngridbz, [eph_fwf, eph_lwf], nfreq_in, temps )
    dirname = string
    interpolation_method = 'spline'
    format = 'text'
    if (associated(input%eph%el_self_energy%output_settings)) then
      if (trim( adjustl( input%eph%el_self_energy%output_settings%directory ) ) /= 'auto') &
        dirname = trim( adjustl( input%eph%el_self_energy%output_settings%directory ) )
      interpolation_method = trim( adjustl( input%eph%el_self_energy%output_settings%interpolation_method ) )
      format = trim( adjustl( input%eph%el_self_energy%output_settings%format ) )
    end if
    if (trim( adjustl( dirname ) ) /= '.') then
      stat = make_directory_pure( dirname, comm=mpiglobal )
      call terminate_if_false( stat==0, '(eph_write_electron_selfenergy_specfun) &
        Failed to create output directory "'//trim( adjustl( dirname ) )//'".' )
    end if
    ! set loop limits
    ik1 = firstofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
    ik2 = lastofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )

    allocate( freqs_out(nfreq_out, fst:lst), selfen_out(nfreq_out, fst:lst) )
    allocate( sfun_in(nfreq_in, fst:lst), sfun_out(nfreq_out, fst:lst) )

    do ik = ik1, ik2
      ! generate self-energy from file
      call eph_else_gen_from_file( selfen_file, ik, freqs_in, temps, el_energy_k, selfen_in, selfen_dw, &
        swidth=input%eph%el_self_energy%swidth )
      if (autofreq) fgrid_out%range = [freqs_in(1), freqs_in(nfreq_in)]

      do itemp = 1, ntemp
        do ist = fst, lst
          ! add FM and DW self-energy
          selfen_in(:, ist, itemp) = selfen_in(:, ist, itemp) + selfen_dw(ist, itemp)
          ! resample self-energy
          fgrid_out%range = fgrid_out%range - el_energy_k(ist)
          call eph_else_resample( nfreq_in, freqs_in-el_energy_k(ist), selfen_in(:, ist, itemp), fgrid_out, freqs_out(:, ist), selfen_out(:, ist), &
            interpolation_method=interpolation_method )
          freqs_out(:, ist) = freqs_out(:, ist) + el_energy_k(ist) 
          fgrid_out%range = fgrid_out%range + el_energy_k(ist)
          ! compute spectral function
          call eph_else_gen_specfun( freqs_in-el_energy_k(ist), selfen_in(:, ist, itemp), sfun_in(:, ist) )
          call eph_else_gen_specfun( freqs_out(:, ist)-el_energy_k(ist), selfen_out(:, ist), sfun_out(:, ist) )
        end do

        ! write self-energy and spectral function to file
        write( fname, '("eph_else+sfun_P",i3.3,"_T",i4.4)' ) ik, nint( temps(itemp) )
        call eph_io_write_else_sfun( reshape( freqs_in, [nfreq_in, 1] ), selfen_in(:, :, itemp), sfun_in, format, trim(dirname)//'/'//trim(fname), el_energy=el_energy_k, selfen_dw=selfen_dw(:, itemp) )
        write( fname, '("eph_else+sfun_rsmpld_P",i3.3,"_T",i4.4)' ) ik, nint( temps(itemp) )
        call eph_io_write_else_sfun( freqs_out, selfen_out, sfun_out, format, trim(dirname)//'/'//trim(fname), el_energy=el_energy_k, selfen_dw=selfen_dw(:, itemp) )
      end do
    end do

    deallocate( selfen_in, selfen_dw, selfen_out, freqs_out, sfun_in, sfun_out )

    ! close binary file
    call selfen_file%close( mpiglobal )
  end subroutine eph_write_electron_selfenergy_specfun
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! ELECTRON QUASI-PARTICLE EQUATION
  !
  !> Solve quasi-particle equation for pre-computed electron self-energy.
  subroutine eph_solve_electron_quasi_particle_equation
    use eph_variables, only: kset => eph_pset ! make target points locally accessible as `kset`
    use eph_inout, only: eph_io_write_energies, eph_io_write_quasi_particle_energies
    use eph_electron_selfenergy

    use block_data_file, only: block_data_file_type
    use mod_kpointset, only: generate_k_vectors
    use mod_eigenvalue_occupancy, only: occmax
    use mod_charge_and_moment, only: chgval
    use mod_mpi_allreduce, only: xmpi_allreduce
    use mod_occupy, only: find_fermi
    use modinput

    integer :: ntemp, fst, lst, ik, ik1, ik2, itemp, ist, stat, i
    real(dp) :: efqp
    character(:), allocatable :: string, format
    character(128) :: fname, dirname
    type(block_data_file_type) :: selfen_file

    real(dp), allocatable :: freqs(:), temps(:), e0(:), selfen_dw(:,:), selfen_fm_hilo(:,:), selfen_dw_hilo(:,:), occ(:,:), &
                             eSP(:,:,:), sDW(:,:,:), sFM_hilo(:,:,:), sDW_hilo(:,:,:)
    complex(dp), allocatable :: selfen_fm(:,:,:), eQP(:,:,:), Z(:,:,:), sFM(:,:,:)
  
    ! set target k-points
    if (associated(input%eph%target)) then
      call eph_var_read_target
    else
      call generate_k_vectors( eph_pset, eph_kset_el%bvec, eph_kset_el%ngridk, eph_kset_el%vkloff, .true., uselibzint=.false. )
    end if
    ! open binary file with electron self-energy
    selfen_file = block_data_file_type( eph_else_filename, [-1], cmplx( 0, 0, dp ) )
    call selfen_file%open( mpiglobal )
    ! read input frequency grid and temperatures from file
    call eph_else_gen_from_file( selfen_file, 1, freqs, temps, e0, selfen_fm, selfen_dw )
    ntemp = size( temps )
    fst = lbound( e0, dim=1 )
    lst = ubound( e0, dim=1 )
    ! set output settings
    string = eph_else_setting_string( input%eph%el_self_energy%ngridbz, [eph_fwf, eph_lwf], size( freqs ), temps )
    dirname = string
    format = 'text'
    if (associated(input%eph%el_self_energy%output_settings)) then
      if (trim( adjustl( input%eph%el_self_energy%output_settings%directory ) ) /= 'auto') &
        dirname = trim( adjustl( input%eph%el_self_energy%output_settings%directory ) )
      format = trim( adjustl( input%eph%el_self_energy%output_settings%format ) )
    end if
    ! set loop limits
    ik1 = firstofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
    ik2 = lastofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )

    allocate( eQP(fst:lst, kset%nkpt, ntemp), Z(fst:lst, kset%nkpt, ntemp), sFM(fst:lst, kset%nkpt, ntemp), source=cmplx( 0, 0, dp ) )
    allocate( eSP(fst:lst, kset%nkpt, ntemp), sDW(fst:lst, kset%nkpt, ntemp), sFM_hilo(fst:lst, kset%nkpt, ntemp), sDW_hilo(fst:lst, kset%nkpt, ntemp), source=0.0_dp )

    do ik = ik1, ik2
      ! generate self-energy from file
      call eph_else_gen_from_file( selfen_file, ik, freqs, temps, e0, selfen_fm, selfen_dw, &
        selfen_fm_hilo=selfen_fm_hilo, &
        selfen_dw_hilo=selfen_dw_hilo, &
        swidth=input%eph%el_self_energy%swidth )
      ! find quasi-particles
      do itemp = 1, ntemp
        do ist = fst, lst
          call eph_else_solve_qp_equation( e0(ist), freqs, selfen_fm(:, ist, itemp)+selfen_dw(ist, itemp), &
            eQP(ist, ik, itemp), sFM(ist, ik, itemp), Z(ist, ik, itemp), &
            input%eph%el_self_energy%eqpsolver )
          eSP(ist, ik, itemp) = e0(ist)
          sDW(ist, ik, itemp) = selfen_dw(ist, itemp)
          sFM(ist, ik, itemp) = sFM(ist, ik, itemp) - sDW(ist, ik, itemp) 
          sFM_hilo(ist, ik, itemp) = selfen_fm_hilo(ist, itemp)
          sDW_hilo(ist, ik, itemp) = selfen_dw_hilo(ist, itemp)
        end do
      end do
    end do
    call xmpi_allreduce( eSP, mpiglobal )
    call xmpi_allreduce( eQP, mpiglobal )
    call xmpi_allreduce( sFM, mpiglobal )
    call xmpi_allreduce( sDW, mpiglobal )
    call xmpi_allreduce( Z, mpiglobal )

    if (mpiglobal%rank == 0) then
      ! write text output
      do itemp = 1, ntemp
        if (eph_pset_path%num_points == 0) then
          write( fname, '("eph_el_qp_energies_T",i4.4)' ) nint( temps(itemp) )
          call eph_io_write_energies( kset%bvec, fst, eQP(:, :, itemp)%re, trim(dirname)//'/'//trim(fname), format, plist=kset%vkl(:, 1:kset%nkpt), &
            xlabel='imaginary part', xvalue=eQP(:, :, itemp)%im )
        else
          write( fname, '("eph_el_qp_disp_T",i4.4)' ) nint( temps(itemp) )
          call eph_io_write_energies( kset%bvec, fst, eQP(:, :, itemp)%re, trim(dirname)//'/'//trim(fname), format, path=eph_pset_path, &
            xlabel='imaginary part', xvalue=eQP(:, :, itemp)%im )
        end if
        write( fname, '("EVALQP_T",i4.4)' ) nint( temps(itemp) )
        call eph_io_write_quasi_particle_energies( kset, fst, eSP(:, :, itemp), eQP(:, :, itemp), sFM(:, :, itemp), sDW(:, :, itemp), sFM_hilo(:, :, itemp), sDW_hilo(:, :, itemp), Z(:, :, itemp), trim(dirname)//'/'//trim(fname), 'text' )
      end do
      ! write binary output
      if (.not. associated(input%eph%target)) then
        allocate( occ(fst:lst, kset%nkpt) )
        do itemp = 1, ntemp
          write( fname, '("EVALQP_T",i4.4,".OUT")' ) nint( temps(itemp) )
          call find_fermi( kset%nkpt, kset%wkpt, lst-fst+1, eQP(:, :, itemp)%re, chgval-occmax*(fst-1), occmax, &
            input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
            efqp, occ )
          call putevalqp( fname, kset, fst, lst, eSP(:, :, itemp), 0.0_dp, eQP(:, :, itemp)%re, efqp )
        end do
      end if
    end if

    ! close binary file
    call selfen_file%close( mpiglobal )
  end subroutine eph_solve_electron_quasi_particle_equation
  !-------------------------------------------------------------------------------- 
  
  !================================================================================ 
  ! ELECTRON EPH COUPLING STRENGTH
  !
  !> Compute electron EPH coupling strength \(\lambda^\pm_{n{\bf k}}\).
  !>
  !> Depending on `parmode`, this routine is either MPI parallelized over target k-points
  !> or over electron bands \(m\) and phonon modes \(\nu\).
  subroutine eph_gen_electron_coupling_strength
    use eph_electrons, only: eph_el_interpolate, eph_el_set_default_frequency_grid
    use eph_phonons, only: eph_ph_energy_q, eph_ph_interpolate, eph_ph_set_default_frequency_grid
    use eph_ephmat, only: eph_ephmat_interpolate
    use eph_inout, only: eph_io_write_coupling_strength, eph_io_write_integrated_coupling_strength
    use eph_eliashberg, only: eph_eliashberg_gen_a2F_tetrahedron

    use mod_kpointset, only: k_set, generate_k_vectors, delete_k_vectors
    use mod_opt_tetra, only: t_set, opt_tetra_destroy, opt_tetra_wgt_delta
    use grid_utils, only: linspace
    use math_utils, only: integrate1d
    use mod_eigenvalue_occupancy, only: occmax
    use exciting_mpi, only: xmpi_allreduce, xmpi_comm_split
    use modinput

    integer :: fmode, lmode, fwf, lwf, ik, ik1, ik2, iq, nfreq, ifreq, ist, n, i
    type(k_set) :: qset
    type(k_set), pointer :: kset
    type(t_set) :: tsetq, tsetk
    type(freq_grid_type) :: fgrid
    type(mpiinfo) :: mpi

    integer, allocatable :: mode_indices(:,:), &  ! phonon mode indices in distributed calculation
                            band_indices(:,:)     ! band indices in distributed calculation
    real(dp), allocatable :: freqs(:), &          ! frequencies
                             vkql(:,:), &         ! shifted BZ integration set k+q
                             el_energy_k(:,:), &  ! electron energies on target set k
                             el_energy_kq(:,:), & ! electron energies on shifted BZ integration set k+q
                             ph_energy_q(:,:), &  ! phonon frequencies on BZ integration set q
                             a2F(:,:,:,:), &      ! Eliashberg spectral function a^2F
                             lambdak(:,:,:), &    ! coupling strength on target set k
                             lambda(:,:), &       ! integrated coupling strength
                             dos(:), &            ! electron density of states
                             wgt(:,:,:)           ! integration weights
    complex(dp), allocatable, target :: el_evec_k(:,:,:), &   ! electron eigenvectors U(k) on target set k
                                        el_evec_kq(:,:,:), &  ! electron eigenvectors U(k+q) on shifted BZ integration set k+q
                                        ph_evec_q(:,:,:), &   ! phonon eigenvectors e(q) on BZ integration set q
                                        ephmat(:,:,:,:,:)     ! EPH matrix elements g(k,q) on BZ integration set
    complex(dp), pointer :: ptr(:,:,:,:)

    ! set k-grid
    if (input%eph%coupling_strength%integrated) then
      allocate( kset )
      call eph_var_init_bz_int( input%eph%coupling_strength%ngridk, input%eph%coupling_strength%vkloff, input%eph%coupling_strength%reducek, kset, tsetk )
    else
      kset => eph_pset
    end if
    ! distribute k-points among processes
    if (input%eph%coupling_strength%parmode == 'points') then
      ik1 = firstofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
      ik2 = lastofset( mpiglobal%rank, kset%nkpt, mpiglobal%procs )
      if (ik1 == 0) then
        ik1 = mod( mpiglobal%rank, kset%nkpt ) + 1
        ik2 = ik1
      end if
      call xmpi_comm_split( mpiglobal, mpi, ik1, mpiglobal%rank )
    else
      mpi = mpiglobal
      ik1 = 1
      ik2 = kset%nkpt
    end if
    ! distribute phonon modes and electron bands among processes
    call patchwork_distribution( eph_nmode, eph_nwf, mpi%procs, mode_indices, band_indices )
    fmode = mode_indices(1, mpi%rank+1) + eph_fmode - 1
    lmode = mode_indices(2, mpi%rank+1) + eph_fmode - 1
    fwf = band_indices(1, mpi%rank+1) + eph_fwf - 1
    lwf = band_indices(2, mpi%rank+1) + eph_fwf - 1
    ! initialize BZ integration points q
    call eph_var_init_bz_int( input%eph%coupling_strength%ngridq, input%eph%coupling_strength%vqloff, input%eph%coupling_strength%reduceq, qset, tsetq )
    ! interpolate electrons on target points k
    call eph_el_interpolate( kset%vkl(:, 1:kset%nkpt), el_energy_k, el_evec_k )
    ! interpolate phonons on BZ integration points q
    call eph_ph_interpolate( qset%vkl(:, 1:qset%nkpt), ph_energy_q, ph_evec_q, mrange=[fmode, lmode] )
    ! set phonon frequency grid
    nfreq = input%eph%coupling_strength%nfreqa2F
    call eph_ph_set_default_frequency_grid( fgrid )
    fgrid%numpoints = nfreq
    freqs = eph_var_gen_frequency_grid( fgrid, reshape( eph_ph_energy_q, [size(eph_ph_energy_q)] ) )

    n = 2; if (input%eph%coupling_strength%quasielastic) n = 1
    allocate( a2F(nfreq, eph_nwf, eph_nwf, n) )
    allocate( lambdak(eph_nwf, kset%nkpt, n), source=0.0_dp )

    ! loop over target points k
    do ik = ik1, ik2
      if (mpi%rank == 0) write( *, '("k ",i4.4,"/",i4.4)' ) ik, kset%nkpt
      ! generate k+q vectors
      vkql = reshape( [(qset%vkl(:, iq)+kset%vkl(:, ik), iq=1, qset%nkpt)], [3, qset%nkpt] )
      ! interpolate electrons on k+q
      call eph_el_interpolate( vkql, el_energy_kq, el_evec_kq, irange=[fwf, lwf], mpicomm=mpi )
      ! interpolate EPH matrix on (k,q)
      ptr(1:size(el_evec_kq, dim=1), 1:size(el_evec_kq, dim=2), 1:1, 1:size(el_evec_kq, dim=3)) => el_evec_kq
      call eph_ephmat_interpolate( kset%vkl(:, ik:ik), qset%vkl(:, 1:qset%nkpt), &
        el_evec_k(:, :, ik:ik), ptr, ph_energy_q, ph_evec_q, &
        ephmat, &
        frange=[fwf, lwf], mrange=[fmode, lmode], mpicomm=mpi )
      ! compute Eliashberg function
      ptr(1:size(ephmat, dim=1), 1:size(ephmat, dim=2), 1:size(ephmat, dim=3), 1:size(ephmat, dim=5)) => ephmat
      if (input%eph%coupling_strength%quasielastic) then
        call eph_eliashberg_gen_a2F_tetrahedron( el_energy_k(:, ik), freqs, el_energy_k(:, ik), el_energy_kq, ph_energy_q, &
          ptr, tsetq, a2F(:, :, :, 1), sig=0 )
      else
        call eph_eliashberg_gen_a2F_tetrahedron( el_energy_k(:, ik), freqs, el_energy_k(:, ik), el_energy_kq, ph_energy_q, &
          ptr, tsetq, a2F(:, :, :, 1), sig=1 )
        call eph_eliashberg_gen_a2F_tetrahedron( el_energy_k(:, ik), freqs, el_energy_k(:, ik), el_energy_kq, ph_energy_q, &
          ptr, tsetq, a2F(:, :, :, 2), sig=-1 )
      end if
      call xmpi_allreduce( a2F, mpi )
      ! compute electron EPH coupling strength
      do ist = 1, eph_nwf
        do i = 1, n
          where (freqs > 2*epsilon(freqs)*maxval(freqs))
            a2F(:, ist, ist, i) = 2 * a2F(:, ist, ist, i) / freqs
          elsewhere
            a2F(:, ist, ist, i) = 0.0_dp
          end where
          lambdak(ist, ik, i) = integrate1d( freqs, a2F(:, ist, ist, i), 'trapez' )
        end do
      end do
    end do
    call xmpi_allreduce( lambdak, mpiglobal )

    ! integrate over k
    if (input%eph%coupling_strength%integrated .and. mpiglobal%rank == 0) then
      ! set electron frequency grid
      call eph_el_set_default_frequency_grid( fgrid )
      fgrid%numpoints = input%eph%coupling_strength%nfreqint
      nfreq = fgrid%numpoints
      freqs = eph_var_gen_frequency_grid( fgrid, reshape( el_energy_k, [kset%nkpt*eph_nwf] ) )
      ! compute electron DOS and integrated coupling strength
      allocate( wgt(eph_nwf, kset%nkpt, nfreq), dos(nfreq), lambda(nfreq, n) )
      call opt_tetra_wgt_delta( tsetk, kset%nkpt, eph_nwf, el_energy_k, nfreq, freqs, wgt )
      do ifreq = 1, nfreq
        dos(ifreq) = sum( wgt(:, :, ifreq) )
        do i = 1, n
          lambda(ifreq, i) = sum( lambdak(:, :, i) * wgt(:, :, ifreq) )
        end do
      end do
      do i = 1, n
        where (dos > 2*epsilon(dos)*maxval(dos))
          lambda(:, i) = lambda(:, i) / dos
        elsewhere
          lambda(:, i) = 0.0_dp
        end where
      end do
      dos = dos * occmax 
    end if

    ! write result to file
    if (input%eph%coupling_strength%integrated) then
      !if (mpiglobal%rank == 0) call eph_io_write_integrated_coupling_strength( freqs, dos, lambda, 'eph_el_integrated_coupling_strength', input%eph%coupling_strength%format )
      if (mpiglobal%rank == 0) call eph_io_write_integrated_coupling_strength( freqs, dos, lambda, 'eph_el_integrated_coupling_strength', 'text' )
      if (mpiglobal%rank == 0) call eph_io_write_integrated_coupling_strength( freqs, dos, lambda, 'eph_el_integrated_coupling_strength', 'json' )
    else
      if (eph_pset_path%num_points == 0) then
        !if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( kset%bvec, eph_fst_span, el_energy_k, lambdak, 'eph_el_coupling_strength', input%eph%coupling_strength%format, plist=kset%vkl(:, 1:kset%nkpt) )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( kset%bvec, eph_fst_span, el_energy_k, lambdak, 'eph_el_coupling_strength', 'text', plist=kset%vkl(:, 1:kset%nkpt) )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( kset%bvec, eph_fst_span, el_energy_k, lambdak, 'eph_el_coupling_strength', 'json', plist=kset%vkl(:, 1:kset%nkpt) )
      else
        !if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( eph_pset%bvec, eph_fst_span, el_energy_k, lambdak, 'eph_el_coupling_strength_disp', input%eph%coupling_strength%format, path=eph_pset_path )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( eph_pset%bvec, eph_fst_span, el_energy_k, lambdak, 'eph_el_coupling_strength_disp', 'text', path=eph_pset_path )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( eph_pset%bvec, eph_fst_span, el_energy_k, lambdak, 'eph_el_coupling_strength_disp', 'json', path=eph_pset_path )
      end if
    end if

    ! clean up
    call delete_k_vectors( qset )
    call opt_tetra_destroy( tsetq )
    if (input%eph%coupling_strength%integrated) then
      call delete_k_vectors( kset )
      call opt_tetra_destroy( tsetk )
    end if
  end subroutine eph_gen_electron_coupling_strength
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! PHONON EPH COUPLING STRENGTH
  !
  !> Compute phonon EPH coupling strength \(\lambda^\pm_{\nu{\bf q}}\).
  !>
  !> Depending on `parmode`, this routine is either MPI parallelized over target q-points
  !> or over electron bands \(m\) and phonon modes \(\nu\).
  subroutine eph_gen_phonon_coupling_strength
    use eph_electrons, only: eph_el_interpolate
    use eph_phonons, only: eph_ph_interpolate, eph_ph_set_default_frequency_grid
    use eph_ephmat, only: eph_ephmat_interpolate
    use eph_inout, only: eph_io_write_coupling_strength, eph_io_write_integrated_coupling_strength
    use eph_eliashberg, only: eph_eliashberg_gen_phonon_coupling_tetrahedron

    use mod_kpointset, only: k_set, generate_k_vectors, delete_k_vectors
    use mod_opt_tetra, only: t_set, opt_tetra_destroy, opt_tetra_wgt_delta, opt_tetra_efermi
    use mod_eigenvalue_occupancy, only: occmax
    use mod_charge_and_moment, only: chgval
    use exciting_mpi, only: xmpi_allreduce, xmpi_comm_split
    use mod_lattice, only: omega
    use unit_conversion, only: bohr_to_m
    use modinput

    integer :: fmode, lmode, fwf, lwf, ik, iq, iq1, iq2, nfreq, ifreq, imode, n, i
    real(dp) :: e0, xchrg, t1
    type(k_set) :: kset
    type(k_set), pointer :: qset
    type(t_set) :: tsetq, tsetk
    type(freq_grid_type) :: fgrid
    type(mpiinfo) :: mpi

    integer, allocatable :: mode_indices(:,:), &  ! phonon mode indices in distributed calculation
                            band_indices(:,:)     ! band indices in distributed calculation
    real(dp), allocatable :: freqs(:), &          ! frequencies
                             vkql(:,:), &         ! shifted BZ integration set k+q
                             el_energy_k(:,:), &  ! electron energies on BZ integration set k
                             el_energy_kq(:,:), & ! electron energies on shifted BZ integration set k+q
                             ph_energy_q(:,:), &  ! phonon frequencies on target set q
                             lambdaq(:,:,:,:), &  ! coupling strength on target set q
                             lambda(:,:), &       ! integrated coupling strength
                             cumlambda(:,:), &    ! cumulative coupling strength
                             dos(:), &            ! phonon density of states
                             wgt(:,:,:), &        ! integration weights
                             f(:), s(:), occ(:,:)
    complex(dp), allocatable, target :: el_evec_k(:,:,:), &   ! electron eigenvectors U(k) on BZ integration set k
                                        el_evec_kq(:,:,:), &  ! electron eigenvectors U(k+q) on shifted BZ integration set k+q
                                        ph_evec_q(:,:,:), &   ! phonon eigenvectors e(q) on target set q
                                        ephmat(:,:,:,:,:)     ! EPH matrix elements g(k,q) on BZ integration set
    complex(dp), pointer :: ptr(:,:,:,:)

    ! set q-grid
    if (input%eph%coupling_strength%integrated) then
      allocate( qset )
      call eph_var_init_bz_int( input%eph%coupling_strength%ngridq, input%eph%coupling_strength%vqloff, input%eph%coupling_strength%reduceq, qset, tsetq )
    else
      qset => eph_pset
    end if
    ! distribute q-points among processes
    if (input%eph%coupling_strength%parmode == 'points') then
      iq1 = firstofset( mpiglobal%rank, qset%nkpt, mpiglobal%procs )
      iq2 = lastofset( mpiglobal%rank, qset%nkpt, mpiglobal%procs )
      if (iq1 == 0) then
        iq1 = mod( mpiglobal%rank, qset%nkpt ) + 1
        iq2 = iq1
      end if
      call xmpi_comm_split( mpiglobal, mpi, iq1, mpiglobal%rank )
    else
      mpi = mpiglobal
      iq1 = 1
      iq2 = qset%nkpt
    end if
    ! distribute phonon modes and electron bands among processes
    call patchwork_distribution( eph_nmode, eph_nwf, mpi%procs, mode_indices, band_indices )
    fmode = mode_indices(1, mpi%rank+1) + eph_fmode - 1
    lmode = mode_indices(2, mpi%rank+1) + eph_fmode - 1
    fwf = band_indices(1, mpi%rank+1) + eph_fwf - 1
    lwf = band_indices(2, mpi%rank+1) + eph_fwf - 1
    ! initialize BZ integration points k
    call eph_var_init_bz_int( input%eph%coupling_strength%ngridk, input%eph%coupling_strength%vkloff, input%eph%coupling_strength%reducek, kset, tsetk )
    ! interpolate phonons on target points q
    call eph_ph_interpolate( qset%vkl(:, 1:qset%nkpt), ph_energy_q, ph_evec_q )
    ! interpolate electrons on BZ integration points k
    call eph_el_interpolate( kset%vkl(:, 1:kset%nkpt), el_energy_k, el_evec_k )
    ! set electron energy to evaluate at
    e0 = input%eph%coupling_strength%energy
    if (input%eph%coupling_strength%doping /= 0.0_dp) then
      xchrg = input%eph%coupling_strength%doping * (bohr_to_m*100)**3 * omega
      if (mpiglobal%rank == 0) print *, 'extra charge:', xchrg
      if (xchrg > epsilon(e0)) then
        t1 = minval( pack( el_energy_k, el_energy_k > 0.0_dp ) )
      else if (xchrg < -epsilon(e0)) then
        t1 = maxval( pack( el_energy_k, el_energy_k < 0.0_dp ) )
      else
        t1 = 0.0_dp
      end if
      if (mpiglobal%rank == 0) print *, 'initial guess energy:', t1
      allocate( occ(eph_fwf:eph_lwf, kset%nkpt) )
      call opt_tetra_efermi( tsetk, (chgval+xchrg)/occmax-eph_fst_span+1, kset%nkpt, eph_nwf, el_energy_k, e0, occ, ef0=t1, df0=0.1_dp )
      deallocate( occ )
      if (mpiglobal%rank == 0) print *, 'evaluation energy:', e0
    end if

    n = 2; if (input%eph%coupling_strength%quasielastic) n = 1
    allocate( lambdaq(1, eph_fmode:eph_lmode, qset%nkpt, n), source=0.0_dp )

    ! loop over target points q
    do iq = iq1, iq2
      if (mpi%rank == 0) write( *, '("q ",i4.4,"/",i4.4)' ) iq, qset%nkpt
      ! generate k+q vectors
      vkql = reshape( [(kset%vkl(:, ik)+qset%vkl(:, iq), ik=1, kset%nkpt)], [3, kset%nkpt] )
      ! interpolate electrons on k+q
      call eph_el_interpolate( vkql, el_energy_kq, el_evec_kq, irange=[fwf, lwf], mpicomm=mpi )
      ! interpolate EPH matrix on (k,q)
      ptr(1:size(el_evec_kq, dim=1), 1:size(el_evec_kq, dim=2), 1:size(el_evec_kq, dim=3), 1:1) => el_evec_kq
      call eph_ephmat_interpolate( kset%vkl(:, 1:kset%nkpt), qset%vkl(:, iq:iq), &
        el_evec_k, ptr, ph_energy_q(fmode:lmode, iq:iq), ph_evec_q(:, fmode:lmode, iq:iq), &
        ephmat, &
        frange=[fwf, lwf], mrange=[fmode, lmode], mpicomm=mpi )
      ! compute coupling strength
      ptr(1:size(ephmat, dim=1), 1:size(ephmat, dim=2), 1:size(ephmat, dim=3), 1:size(ephmat, dim=4)) => ephmat
      if (input%eph%coupling_strength%quasielastic) then
        call eph_eliashberg_gen_phonon_coupling_tetrahedron( [e0], el_energy_k, el_energy_kq, ph_energy_q(fmode:lmode, iq), &
          ptr, tsetk, lambdaq(:, fmode:lmode, iq, 1), sig=0 )
      else
        call eph_eliashberg_gen_phonon_coupling_tetrahedron( [e0], el_energy_k, el_energy_kq, ph_energy_q(fmode:lmode, iq), &
          ptr, tsetk, lambdaq(:, fmode:lmode, iq, 1), sig=1 )
        call eph_eliashberg_gen_phonon_coupling_tetrahedron( [e0], el_energy_k, el_energy_kq, ph_energy_q(fmode:lmode, iq), &
          ptr, tsetk, lambdaq(:, fmode:lmode, iq, 2), sig=-1 )
      end if
    end do
    call xmpi_allreduce( lambdaq, mpiglobal )

    ! integrate over q
    if (input%eph%coupling_strength%integrated .and. mpiglobal%rank == 0) then
      ! set phonon frequency grid
      call eph_ph_set_default_frequency_grid( fgrid )
      fgrid%numpoints = input%eph%coupling_strength%nfreqint
      nfreq = fgrid%numpoints
      freqs = eph_var_gen_frequency_grid( fgrid, reshape( ph_energy_q, [qset%nkpt*eph_nmode] ) )
      ! compute phonon DOS and integrated coupling strength
      allocate( wgt(eph_nmode, qset%nkpt, nfreq), dos(nfreq), lambda(nfreq, n), cumlambda(nfreq, n) )
      allocate( f(nfreq) )
      call opt_tetra_wgt_delta( tsetq, qset%nkpt, eph_nmode, ph_energy_q, nfreq, freqs, wgt )
      do ifreq = 1, nfreq
        dos(ifreq) = sum( wgt(:, :, ifreq) )
        do i = 1, n
          lambda(ifreq, i) = sum( lambdaq(1, :, :, i) * ph_energy_q * wgt(:, :, ifreq) )
        end do
      end do
      ! compute cumulated coupling strength
      do i = 1, n
        where (freqs > 2*epsilon(freqs)*maxval(freqs))
          f = 2 * lambda(:, i) / freqs
        elsewhere
          f = 0.0_dp
        end where
        s = (freqs(2:) - freqs(:nfreq-1)) * (f(2:) + f(:nfreq-1)) / 2
        cumlambda(1, i) = 0.0_dp
        cumlambda(2:, i) = [(sum(s(:ifreq)), ifreq=1, size(s))]
      end do
    end if

    ! write result to file
    if (input%eph%coupling_strength%integrated) then
      !if (mpiglobal%rank == 0) call eph_io_write_integrated_coupling_strength( freqs, dos, lambda, 'eph_ph_integrated_coupling_strength', input%eph%coupling_strength%format, &
      !  cumulative_lambda=cumlambda )
      if (mpiglobal%rank == 0) call eph_io_write_integrated_coupling_strength( freqs, dos, lambda, 'eph_ph_integrated_coupling_strength', 'text', &
        cumulative_lambda=cumlambda )
      if (mpiglobal%rank == 0) call eph_io_write_integrated_coupling_strength( freqs, dos, lambda, 'eph_ph_integrated_coupling_strength', 'json', &
        cumulative_lambda=cumlambda )
    else
      if (eph_pset_path%num_points == 0) then
        !if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( qset%bvec, eph_fmode, ph_energy_q, lambdaq(1,:,:,:), 'eph_ph_coupling_strength', input%eph%coupling_strength%format, plist=qset%vkl(:, 1:qset%nkpt) )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( qset%bvec, eph_fmode, ph_energy_q, lambdaq(1,:,:,:), 'eph_ph_coupling_strength', 'text', plist=qset%vkl(:, 1:qset%nkpt) )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( qset%bvec, eph_fmode, ph_energy_q, lambdaq(1,:,:,:), 'eph_ph_coupling_strength', 'json', plist=qset%vkl(:, 1:qset%nkpt) )
      else
        !if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( eph_pset%bvec, eph_fmode, ph_energy_q, lambdaq(1,:,:,:), 'eph_ph_coupling_strength_disp', input%eph%coupling_strength%format, path=eph_pset_path )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( eph_pset%bvec, eph_fmode, ph_energy_q, lambdaq(1,:,:,:), 'eph_ph_coupling_strength_disp', 'text', path=eph_pset_path )
        if (mpiglobal%rank == 0) call eph_io_write_coupling_strength( eph_pset%bvec, eph_fmode, ph_energy_q, lambdaq(1,:,:,:), 'eph_ph_coupling_strength_disp', 'json', path=eph_pset_path )
      end if
    end if

    ! clean up
    call delete_k_vectors( kset )
    call opt_tetra_destroy( tsetk )
    if (input%eph%coupling_strength%integrated) then
      call delete_k_vectors( qset )
      call opt_tetra_destroy( tsetq )
    end if
  end subroutine eph_gen_phonon_coupling_strength
  !-------------------------------------------------------------------------------- 
end module eph
