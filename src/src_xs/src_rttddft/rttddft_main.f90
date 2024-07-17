! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created Apr 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Cleaned May 2023 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module is the kernel of a RT-TDDFT calculation.
!> It contains the subroutine `coordinate_rttddft_calculation`, which manages
!> a RT-TDDFT calculation. Also here are implemented the following subroutines:
!> uppot, uprho.
module rttddft_main
  use asserts, only: assert
  use m_getunit, only: getunit
  use MD, only: force, MD_input_keys
  use MD_io, only: MD_out
  use mod_atoms, only: natmtot, natoms, nspecies, atposc, idxas
  use mod_charge_and_moment, only: chgval
  use mod_eigenvalue_occupancy, only: occsv, nstfv
  use mod_eigensystem, only: nmatmax, nmat
  use mod_kpoint, only: nkpt
  use mod_lattice, only: omega
  use mod_misc, only: filext
  use mod_mpi_env, only: mpiinfo
  use modinput, only: input, input_type
  use modmpi, only: rank, procs, mpi_env_k, mpiglobal, distribute_loop, barrier, terminate_if_false
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_CurrentDensity, only: Current_Density_Paramagnetic_Compoment
  use rttddft_Density, only: UpdateDensity
  use rttddft_Energy, only: TotalEnergy, obtain_energy_rttddft
  use rttddft_GlobalVariables
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_init, only: initialize_rttddft
  use rttddft_input, only: rttddft_input_keys
  use rttddft_io, only: open_files_jpa, close_files_jpa, write_jpa, &
    open_file_timing, close_file_timing, write_timing, &
    open_file_nexc, close_file_nexc, write_nexc, &
    open_file_etot, close_file_etot, write_total_energy, &
    open_file_info, close_file_info, write_file_info, write_file_info_header, &
    write_wavefunction
  use rttddft_MD, only: force_rttdft, move_ions, &
    MD_allocate_global_arrays => allocate_global_arrays, &
    MD_deallocate_global_arrays => deallocate_global_arrays, &
    MD_evaluate_charge_val => evaluate_charge_val
  use rttddft_NumberExcitations, only: Obtain_number_excitations
  use rttddft_screenshot, only: screenshot
  use rttddft_timings, only: Timing_RTTDDFT_and_MD, Timing_RTTDDFT, &
    Timing_RTTDDFT_density, Timing_RTTDDFT_potential, Print_Timings, timesec_RTTDDFT
  use rttddft_VectorPotential, only: solver_types, Calculate_Vector_Potential, Evolve_A_ind => Solve_ODE_Vector_Potential
  use rttddft_Wavefunction, only: UpdateWavefunction, Update_basis_derivative, SE, EH, propagator_types
  
  implicit none

  private
  
  public  :: coordinate_rttddft_calculation

contains

  !> This subroutine manages a RT-TDDFT calculation.
  !> <ol>
  !> <li> Run a single-shot groundstate calculation using the already converged
  !> density and potential. </li>
  !> <li> Obtain the KS wavefunctions, the density and the hamiltonian at 
  !> \( t=0 \). </li>
  !> <li> Evolve the wavefunctions, the density and the hamiltonian using the 
  !> desired time step. </li>
  !> </ol>
  subroutine coordinate_rttddft_calculation()
    ! counter for the number of iterations of real-time
    integer                 :: it
    ! indexes of the first and the last k-points
    integer                 :: first_kpt, last_kpt
    ! Number of time steps \( \Delta t \) required to reach `tend`
    integer                 :: n_steps

    integer                 :: iprint
    integer                 :: is, ia, ias, timeStepMultiplier, l_rad_step
    logical                 :: predCorrReachedMaxSteps

    character(len=100)      :: string

    real(dp), allocatable   :: atom_positions(:, :) ! in cartesian coordinates x, y, z
    real(dp), allocatable   :: atom_velocities(:, :) ! in cartesian coordinates x, y, z
    type(force)             :: forces
    type(MD_input_keys)     :: molecular_dynamics
    type(rttddft_input_keys):: rt

    ! Current time \( t \) for the time evolution carried out in RT-TDDFT
    real(dp)                :: time

    real(dp),allocatable    :: nex(:), ngs(:), nt(:)
    real(dp)                :: aindsave(3),pvecsave(3)
    real(dp)                :: jindsave(3),aextsave(3),atotsave(3)
    real(dp)                :: electric_field(3)
    real(dp)                :: timei, timef, timeiter
    real(dp)                :: tol
    real(dp), parameter     :: tol_default = 1e-10_dp
    type(MD_out)            :: MD_outputs

    ! Variables to store data and print
    real(dp),allocatable    :: timestore(:),aindstore(:,:),atotstore(:,:)
    real(dp),allocatable    :: jindstore(:,:), pvecstore(:,:)
    real(dp),allocatable    :: atposcstore(:,:,:), velstore(:,:,:)
    type(force),allocatable :: forces_store(:)
    logical,allocatable     :: printforces(:)
    logical, allocatable    :: screenshot_was_taken(:)

    type(TotalEnergy), allocatable  :: etotstore(:)
    type(Print_Timings)             :: printTimings
    type(Timing_RTTDDFT_and_MD)     :: timing
    type(Timing_RTTDDFT_and_MD), allocatable :: timing_store(:)


    call timesec( timei )

    ! Sanity check
    call sanity_checks( input )

    ! Interface with input parameters
    tol = tol_default
    if( associated(input%groundstate%solver) ) tol = input%groundstate%solver%evaltol
    call rt%parse_input( input%xs%realTimeTDDFT, tol )
    call molecular_dynamics%parse_input()
    time = 0._dp
    n_steps = int( rt%t_end / rt%propagator%time_step )
    l_rad_step = input%groundstate%lradstep
    
    ! we only perform MD in RT-TDDFT if the type is Ehrenfest
    if( molecular_dynamics%on ) molecular_dynamics%on = ( trim(molecular_dynamics%MD_type) == 'Ehrenfest' )
    
    ! Outputs general info to RTTDDFT_INFO.OUT 
    if( rank == 0 ) then
      call open_file_info
      call write_file_info_header
    end if

    call initialize_rttddft( rt%pmat, rt%predictor_corrector%on, molecular_dynamics )
    
    if( molecular_dynamics%on ) call init_MD( time, rt%propagator%time_step, timeStepMultiplier, molecular_dynamics, &
      MD_outputs, atom_positions, atom_velocities, electric_field, forces )
    
    call printTimings%set( rt%timings_general, rt%timings_detailed )

    ! Allocate variables to be stored and printed only after rt_input%n_print steps
    allocate(timestore(rt%n_print), aindstore(3,rt%n_print), atotstore(3,rt%n_print))
    allocate(jindstore(3,rt%n_print), pvecstore(3,rt%n_print))
    if( printTimings%general() ) then
      allocate( timing_store(rt%n_print) )
      allocate( screenshot_was_taken(rt%n_print), source=.False. )
    end if
    if( rt%calculate_total_energy ) allocate(etotstore(rt%n_print))
    if( rt%calculate_n_exc ) allocate(nex(rt%n_print),ngs(rt%n_print),nt(rt%n_print))
    if( molecular_dynamics%on ) then
      allocate( printforces(rt%n_print), atposcstore(3,natmtot,rt%n_print), velstore(3,natmtot,rt%n_print))
      if ( molecular_dynamics%print_all_force_components ) then
        allocate( forces_store(rt%n_print) )
        do is = 1, rt%n_print
          call forces_store(is)%allocate_arrays( natmtot )
        end do
      end if
    end if

    if( rank == 0 ) then
      call open_files_jpa
      call write_jpa( time, aind, atot, label='avec' )
      call write_jpa( time, pvec, label='pvec' )
      call write_jpa( time, jind, label='jind' )
    end if

    ! Initialize integers that contain the first and last k-point
    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    ! Total energy
    if ( rt%calculate_total_energy ) then
      call potcoul
      call potxc
      call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_gnd, mpi_env_k, etotstore(1) )
      ! Trick: we need an array to call the subroutine print_total_energy
      timestore(1) = time
      if( rank == 0 ) then
        call open_file_etot
        call write_total_energy( .True., 1, timestore(1), etotstore(1) )
      end if
    end if

    ! Number of excitations
    if (rt%calculate_n_exc) then
      call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
        & evecfv_time, overlap, mpi_env_k, nex(1), ngs(1), nt(1) )
      ! Trick: we need an array to call the subroutine print_nexc
      timestore(1) = time
      if( rank == 0 ) then
        call open_file_nexc
        call write_nexc( .True., 1, timestore(1), nex(1), ngs(1), nt(1) )
      end if
    end if

    if ( rt%screenshots%on ) call screenshot( 0, first_kpt, last_kpt, overlap, evecfv_gnd, &
        & evecfv_time, ham_time )

    if( printTimings%general() ) then
      call timesec( timef )
      if( rank == 0 ) then 
        call open_file_timing
        call write_timing( timef-timei ) !write time for initialization
      end if
    end if

    iprint = 1
    timeiter = timef
    ! This is the most important loop (performed for each time step \(\Delta t\)
    do it = 1, n_steps
      ! Variable to store the timing of each iteration
      timei = timeiter

      ! The "real time" t of our evolution
      time = time + rt%propagator%time_step

      ! WAVEFUNCTION
      if ( rt%predictor_corrector%on ) evecfv_save(:,:,:) = evecfv_time(:,:,:)
      if ( molecular_dynamics%on .and. molecular_dynamics%basis_derivative ) then
        call UpdateWavefunction( first_kpt, rt%propagator, .False., &
        ham_time(:, :, first_kpt : last_kpt), &
        ham_past(:, :, first_kpt : last_kpt), &
        evecfv_time(:, :, first_kpt : last_kpt), &
        overlap(:, :, first_kpt : last_kpt), nmat(1, first_kpt : last_kpt ), &
        atom_velocities )
      else 
        call UpdateWavefunction( first_kpt, rt%propagator, .False., &
        ham_time(:, :, first_kpt : last_kpt), &
        ham_past(:, :, first_kpt : last_kpt), &
        evecfv_time(:, :, first_kpt : last_kpt), &
        overlap(:, :, first_kpt : last_kpt), nmat(1, first_kpt : last_kpt ) )
      end if
      if ( printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%wavefunction )

      ! Update the paramagnetic component of the induced current density
      jparanext = Current_Density_Paramagnetic_Compoment( evecfv_time, pmat, occsv(:, first_kpt:last_kpt), &
        [(1._dp/nkpt, is = first_kpt, last_kpt)], mpi_env_k )
      if ( rt%subtract_J0 ) jparanext(:) = jparanext(:)-jparaspurious(:)
      if ( printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%current_density )

      ! DENSITY
      call uprho( it, rt%propagator%normalize_WF, l_rad_step, printTimings, timing%t_RTTDDFT%dens )

      ! KS-POTENTIAL
      call uppot( printTimings, timing%t_RTTDDFT%pot )

      ! VECTOR POTENTIAL
      if( printTimings%general() ) call timesec( timei )
      ! Check if we need to save aind, pvec, atot and aext
      if( rt%is_field_type_external() .and. rt%predictor_corrector%on .and. ( .not. rt%is_solver_euler() ) ) then
        aindsave(:) = aind(:)
        pvecsave(:) = pvec(:)
        aextsave(:) = aext(:)
      end if
      atotsave = atot
      if( rt%is_field_type_total() ) then
        call update_vector_potential( time, rt%propagator%time_step, rt%vector_potential_solver, atot )
        if( molecular_dynamics%on ) then
          call Calculate_Vector_Potential( time+rt%propagator%time_step, aindsave ) ! trick: aindsave is an auxiliary variable
          electric_field = obtain_electric_field( 2*rt%propagator%time_step, aindsave, atotsave )
        end if
      else 
        call update_vector_potential( time, rt%propagator%time_step, rt%vector_potential_solver, atot, aind, aext )
        if( molecular_dynamics%on ) electric_field = obtain_electric_field( rt%propagator%time_step, atot, atotsave )
      end if
      if( printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%vector_potential )

      ! INDUCED CURRENT
      ! Update the diamagnetic component of the induced current density
      jdia(:) = - atot(:) * chgval / c / omega
      ! Update the paramagnetic component of the induced current density
      jparaold(:) = jpara(:)
      jpara(:) = jparanext(:)
      ! Update the total induced current
      if ( rt%predictor_corrector%on .and. ( .not. rt%is_solver_euler() ) ) then
        jindsave(:) = jind(:)
      end if
      jind(:) = jpara(:) + jdia(:)

      ! HAMILTONIAN
      call UpdateHam( predcorr=.False., calculateOverlap=.False., forcePmatHermitian=rt%pmat%force_pmat_hermitian, &
        printTimings=printTimings, t_ham=timing%t_RTTDDFT%ham, t_MD=timing%t_Ehrenfest, &
        update_mathcalH=.False., update_mathcalB=.False., update_pmat=.False. )

      ! Remark: it makes no sense to employ the predictor-corrector method with SE or EH!
      if ( rt%predictor_corrector%on .and. (rt%propagator%name /= SE) .and. (rt%propagator%name /= EH) ) then
        if ( printTimings%general() ) call timesec( timei )
        call loopPredictorCorrector( it, time, rt, l_rad_step, rt%is_field_type_external() .and. ( .not. rt%is_solver_euler() ), &
          first_kpt, last_kpt, aindsave, atotsave, aextsave, pvecsave, jindsave, &
          jparaold, mpi_env_k, predCorrReachedMaxSteps )
        
        if ( predCorrReachedMaxSteps .and. rank == 0 ) &
          write(*,*) 'Problems with convergence (PredCorr), time: ', time
        if ( molecular_dynamics%on .and. rt%is_field_type_external() ) &
          electric_field(:) = (-1.0_dp/c/rt%propagator%time_step)*(atot(:)-atotsave(:))
        if ( printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%pred_corr )
      end if !predictor-corrector

      ! Obtain the total energy, if requested
      if( rt%calculate_total_energy ) then
        if ( printTimings%detailed() ) call timesec( timei )
        call obtain_energy_rttddft( first_kpt, last_kpt, ham_time, evecfv_time, mpi_env_k, etotstore(iprint) )
        if ( printTimings%detailed() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%energy )
      end if

      ! Obtain the number of excited electrons, if requested
      if( rt%calculate_n_exc ) then
        if ( printTimings%detailed() ) call timesec( timei )
        call Obtain_number_excitations( first_kpt, last_kpt, evecfv_gnd, &
          & evecfv_time, overlap, mpi_env_k, nex(iprint), ngs(iprint), nt(iprint))
        if( printTimings%detailed() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%n_exc )
      end if

      if ( molecular_dynamics%on ) then
        if ( mod( it, timeStepMultiplier ) == 0 ) then
          if ( printTimings%general() ) then 
            call timesec( timei )
            timing%t_Ehrenfest%MD_was_carried_out = .True.
          end if
          call forces%save_total_force()
          call force_rttdft( forces, electric_field, molecular_dynamics, printTimings, timing%t_Ehrenfest )
          call move_ions( forces%total, forces%total_save, molecular_dynamics%time_step, &
            atom_velocities, printTimings, timing%t_Ehrenfest )
          printforces(iprint) = .True.
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia,is)
              atposcstore(1:3, ias, iprint) = atposc(1:3, ia, is)
            end do
          end do
          velstore(:,:,iprint) = atom_velocities(:,:)
          forces_store(iprint) = forces
          ! Update Hamiltonian with the new basis
          if( molecular_dynamics%update_overlap .or. allocated(mathcalH) .or. &
            & allocated(mathcalB) .or. molecular_dynamics%update_pmat ) then
            call UpdateHam( predcorr=.False., forcePmatHermitian=rt%pmat%force_pmat_hermitian, &
              & calculateOverlap=molecular_dynamics%update_overlap, &
              & printTimings=printTimings, t_ham=timing%t_RTTDDFT%ham, t_MD=timing%t_Ehrenfest, &
              & update_mathcalH=allocated(mathcalH), &
              & update_mathcalB=allocated(mathcalB), &
              & update_pmat=molecular_dynamics%update_pmat )
          end if
          if( printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_Ehrenfest%t_MD_step )
        else ! if ( mod( it, timeStepMultiplier ) == 0 )
          printforces(iprint) = .False.
          if ( printTimings%general() ) timing%t_Ehrenfest%MD_was_carried_out = .False.
        end if ! if ( mod( it, timeStepMultiplier ) == 0 )
      end if ! if ( molecular_dynamics%on ) then

      ! Check if a screenshot has been requested
      if ( rt%screenshots%on ) then
        if ( mod( it, rt%screenshots%n_steps ) == 0 ) then
          if( printTimings%general() ) screenshot_was_taken(iprint) = .True.
          if( printTimings%general() ) call timesec(timei)
          call screenshot( it, first_kpt, last_kpt, overlap, evecfv_gnd, &
            & evecfv_time, ham_time )
          if( printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%screenshot )
        else 
          if( printTimings%general() ) screenshot_was_taken(iprint) = .False.
        end if
      end if

      ! Store relevant information from this iteration
      timestore(iprint) = time
      aindstore(:, iprint) = aind(:)
      atotstore(:, iprint) = atot(:)
      pvecstore(:, iprint) = pvec(:)
      jindstore(:, iprint) = jind(:)
      if( printTimings%general() ) timing_store(iprint) = timing

      ! Print relevant information, every 'rt_input%n_print' steps
      if ( iprint == rt%n_print ) then
        ! Update the counter
        iprint = 1
        if( rank == 0 ) then
          call write_jpa( timestore, aindstore, atotstore, label='avec' )
          call write_jpa( timestore, pvecstore, label='pvec' )
          call write_jpa( timestore, jindstore, label='jind' )
          if ( rt%calculate_total_energy ) call write_total_energy( .False., rt%n_print, &
            timestore(:), etotstore(:) )
          if ( rt%calculate_n_exc ) call write_nexc( .False., rt%n_print, timestore(:), &
            nex(:), ngs(:), nt(:) )

          ! Print forces - if this has been requested
          if( molecular_dynamics%on ) then
            do iprint = 1, rt%n_print
              if( printforces(iprint) ) call write_MD_outputs( timestore(iprint), &
                atposcstore(:, :, iprint), velstore(:,:,iprint), forces_store(iprint), &
                molecular_dynamics%print_all_force_components, MD_outputs )
            end do
          end if ! if( molecular_dynamics%on )
        end if

        ! Update the counter
        iprint = 1
        if( printTimings%general() ) then
          call timesec_RTTDDFT( timeiter, timing_store(rt%n_print)%t_iteration )
          call write_timing( it, printTimings%detailed(), rt%calculate_total_energy, &
            rt%calculate_n_exc, rt%predictor_corrector%on, timing_store, screenshot_was_taken, molecular_dynamics%on )
        end if
      else ! if ( iprint .eq. rt_input%n_print ) then
        if( printTimings%general() ) call timesec_RTTDDFT( timeiter, timing_store(iprint)%t_iteration )
        iprint = iprint + 1
      end if ! if ( iprint == rt_input%n_print ) 
      ! Make all the processes wait here: the master alone has been writing the files above
      call barrier( mpi_env_k )
    end do ! do it = 1, nsteps

    if ( rank == 0 ) then
      call close_files_jpa
      call write_file_info( 'Real-time TDDFT calculation finished' )
      call close_file_info
      if( rt%calculate_total_energy ) call close_file_etot
      if( rt%calculate_n_exc ) call close_file_nexc
      if( printTimings%general() ) call close_file_timing
      if( molecular_dynamics%on ) call MD_outputs%close_files()
    end if

    ! write wavefunction, and potential and density with _RTTDDFT.OUT as suffix
    string = filext
    filext = '_RTTDDFT'//trim(filext)
    call write_wavefunction( first_kpt, evecfv_time )
    if ( rank == 0 ) call writestate
    filext = string

    call deallocate_global_arrays( rt%predictor_corrector%on, molecular_dynamics%on )
    
  end subroutine coordinate_rttddft_calculation

  !> This is just an interface to call the subroutine `[[UpdateDensity]]`, which
  !> updates the charge density
  subroutine uprho( iteration_counter, normalize, l_rad_step, printTimings, timing_dens )
    !> Tells how many time steps have already been executed
    integer, intent(in) :: iteration_counter
    !> If `.true.`, normalize the charge density
    logical, intent(in)             :: normalize
    !> radial step length
    integer(i32), intent(in)        :: l_rad_step
    !> Object that packs information about printing of timings [[Print_Timings]]
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(inout) :: timing_dens

    if( present(printTimings) ) then
      call assert( present(timing_dens), 'timing_dens must be also present when printTimings is' )
      call UpdateDensity( iteration_counter, normalize, l_rad_step, printTimings, timing_dens )
    else
      call UpdateDensity( iteration_counter, normalize, l_rad_step )
    end if
  end subroutine uprho


  !> This is just an interface to call the subroutines that updates the KS potential
  subroutine uppot( printTimings, t_pot )
    !> Object that packs information about printing of timings [[Print_Timings]]
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the KS potential
    type(Timing_RTTDDFT_potential), optional, intent(inout) :: t_pot

    logical :: timings_general, timings_detailed
    real(dp) :: ti, tstart

    timings_general = .false.
    timings_detailed = .false.
    if( present(printTimings) ) then 
      call assert( present(t_pot), 't_pot must be present if printTimings is')
      call printTimings%get( timings_general, timings_detailed )
    end if 
    if( timings_general ) then
      call timesec( ti )
      tstart = ti
    end if

    ! Compute the effective potential (with the updated density)
    call poteff
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_pot%poteff )
    ! Fourier transform effective potential to G-space
    call genveffig
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_pot%genveffig )
    call genmeffig
    if( timings_general ) then
      call timesec_RTTDDFT( tstart, t_pot%total )
      if( timings_detailed ) t_pot%genmeffig = tstart-ti
    end if
  end subroutine uppot

  !> (private subroutine) Check if variables given in the input file make sense
  subroutine sanity_checks( inp  )
    !> type with the variables given in the input file
    type(input_type):: inp
    
    call terminate_if_false( .not. inp%groundstate%solver%packedmatrixstorage, &
      & 'Error: RT-TDDFT does not work with matrices stored in a packed form.' )

    ! Consistency check: check if no spin polarized calculations are requested.
    call terminate_if_false( .not. inp%groundstate%tevecsv, &
      & 'Error: only spin unpolarised calculations are possible with RT-TDDFT now.' )

    ! Consistency check: laser has been defined?
    call terminate_if_false( associated(inp%xs%realTimeTDDFT%laser), &
      & 'Element <laser> in <realTimeTDDFT> not found')

    ! Consistency check
    call terminate_if_false( associated(inp%xs%realTimeTDDFT%pmat), &
      & 'Element <pmat> in <realTimeTDDFT> not found' )

    ! Consistency check: MD and predictor corrector?
    call terminate_if_false( .not. ( associated(inp%xs%realTimeTDDFT%predictorCorrector) .and. associated(inp%MD) ), &
      & 'It is currently not possible to use the predictor corrector method together with molecular dynamics')

  end subroutine

  subroutine update_vector_potential( t, dt, method, A_tot, A_ind, A_ext )
    !> time \(t\)
    real(dp), intent(in)              :: t
    !> time step \( \Delta t\)
    real(dp), intent(in)              :: dt
    !> Method used to update the vector potential
    integer(kind(solver_types)), intent(in) :: method
    !> total vector potential: \(A_{tot} = A_{ind} + A_{ext}\)
    real(dp), intent(inout)           :: A_tot(3)
    !> induced vector potential
    real(dp), intent(inout), optional :: A_ind(3)
    !> external vector potential
    real(dp), intent(inout), optional :: A_ext(3)

    logical :: all_fields_present

    all_fields_present = present(A_ind)
    if( all_fields_present ) call assert( present(A_ext), 'If A_ind is present, then A_ext must also be' )

    !TODO(Ronaldo): Refactor `Evolve_A_ind` to avoid globals
    call Evolve_A_ind( t, dt, method, .not. all_fields_present )
    if( all_fields_present ) then
      call Calculate_Vector_Potential( t, A_ext(:) )
      A_tot(:) = A_ind(:) + A_ext(:)
    else
      call Calculate_Vector_Potential( t, A_tot(:) )
    end if
  end subroutine

  !> Obtain the electric field as the time derivative of the vector potential
  pure function obtain_electric_field( dt, A_tot, A_tot_previous ) result( E_field )
    !> time step
    real(dp), intent(in)  :: dt
    !> vector potential at time `t`
    real(dp), intent(in)  :: A_tot(3)
    !> vector potential at time `t-dt`
    real(dp), intent(in)  :: A_tot_previous(3)
    real(dp) :: E_field(3)
    E_field = ( -1._dp / c / dt ) * ( A_tot - A_tot_previous )
  end function

  subroutine loopPredictorCorrector( it, time, rt, l_rad_step, evolveA, first_kpt, last_kpt, &
    aindsave, atotsave, aextsave, pvecsave, jindsave, jparasave, mpi_env, maxStepsReached )
    !> current iteration number in the RT-TDDFT loop
    integer(i32), intent(in)       :: it
    !> time \( t \)
    real(dp), intent(in)  :: time
    !> Type that encapsulates the parameters defined in the input file
    type(rttddft_input_keys), intent(in) :: rt
    !> radial step length
    integer(i32), intent(in)        :: l_rad_step
    !> If `.True.`, the vector potential is evolved in each step
    logical, intent(in)            :: evolveA
    !> index of the first `k-point` to be considered in the sum
    integer(i32),intent(in)        :: first_kpt
    !> index of the last `k-point` considered
    integer(i32),intent(in)        :: last_kpt
    !> Backup of `aind`
    real(dp), intent(in)           :: aindsave(3)
    !> Backup of `atot`
    real(dp), intent(in)           :: atotsave(3)
    !> Backup of `aext`
    real(dp), intent(in)           :: aextsave(3)
    !> Backup of `pvec`
    real(dp), intent(in)           :: pvecsave(3)
    !> Backup of `jind`
    real(dp), intent(in)           :: jindsave(3)
    !> Backup of `jpara`
    real(dp), intent(in)           :: jparasave(3)
    !> MPI environment
    type(mpiinfo), intent(in)      :: mpi_env
    !> When `.True.`, it informs that the maximum steps have been reached
    logical, intent(out)           :: maxStepsReached

    integer(i32) :: i, ik
    real(dp)     :: err

    do i = 1, rt%predictor_corrector%max_steps
      ! WAVEFUNCTION
      evecfv_time(:,:,:) = evecfv_save(:,:,:)
      call UpdateWavefunction( first_kpt, rt%propagator, .True., &
      ham_time(:, :, first_kpt : last_kpt), &
      ham_past(:, :, first_kpt : last_kpt), &
      evecfv_time(:, :, first_kpt : last_kpt), &
      overlap(:, :, first_kpt : last_kpt), nmat(1, first_kpt : last_kpt ) )

      ! Update the paramagnetic component of the induced current density
      jparanext = Current_Density_Paramagnetic_Compoment( evecfv_time, pmat, occsv(:, first_kpt:last_kpt), &
        [(1._dp/nkpt, ik = first_kpt, last_kpt)], mpi_env )
      if ( rt%subtract_J0 ) jparanext(:) = jparanext(:)-jparaspurious(:)

      ! DENSITY
      call uprho( it, rt%propagator%normalize_WF, l_rad_step )

      ! KS-POTENTIAL
      call uppot

      ! VECTOR POTENTIAL
      ! Update the induced part of the vector potential
      if( evolveA ) then
        jpara(:) = jparasave(:) !attention: jparaold saves the value of jpara(t-deltat)
        aind(:) = aindsave(:)
        atot(:) = atotsave(:)
        aext(:) = aextsave(:)
        pvec(:) = pvecsave(:)
        jind(:) = jindsave(:)
        call update_vector_potential( time, rt%propagator%time_step, rt%vector_potential_solver, atot, aind, aext )
        call Evolve_A_ind( time, rt%propagator%time_step, rt%vector_potential_solver, .False. )
        call Calculate_Vector_Potential( time, aext(:) )
        ! Update the (total) vector potential
        atot(:) = aind(:) + aext(:)
      end if

      ! INDUCED CURRENT
      ! Update the paramagnetic component of the induced current density
      jdia(:) = -atot(:)*chgval/c/omega
      ! Update the paramagnetic component of the induced current density
      jpara(:) = jparanext(:)
      jind(:) = jpara(:)+jdia(:)

      ! HAMILTONIAN
      ham_predcorr(:,:,:) = ham_time(:,:,:)
      call UpdateHam( predcorr=.True., calculateOverlap=.False., forcePmatHermitian=rt%pmat%force_pmat_hermitian )

      ! Check the difference between the two hamiltonians
      err = maxval(abs(ham_predcorr(:,:,:)-ham_time(:,:,:)))
      if ( err <= rt%predictor_corrector%tol ) exit

    end do
    maxStepsReached = (i>rt%predictor_corrector%max_steps)
  end subroutine 

  !> Subroutine to initialize all MD related variables
  subroutine init_MD( t_0, timeStepRTTDDFT, timeStepMultiplier, molecular_dynamics, &
      MD_outputs, atom_positions, atom_velocities, e_field, forces )
    !> Initial time \( t_0 \)
    real(dp), intent(in)               :: t_0
    !> Time step used in the real-time TDDFT calculation
    real(dp), intent(in)               :: timeStepRTTDDFT
    !> Integer ratio between the time step used in MD and `timeStepRTTDDFT`
    integer(i32), intent(out)          :: timeStepMultiplier
    !> variable with interfaces to elements defined in the input file
    type(MD_input_keys), intent(inout) :: molecular_dynamics
    !> variable with interfaces to MD outputs
    type(MD_out), intent(inout)        :: MD_outputs
    !> positions of all atoms in cartesian coordinates
    real(dp), allocatable, intent(out) :: atom_positions(:, :) 
    !> velocities of all atoms in cartesian coordinates
    real(dp), allocatable, intent(out) :: atom_velocities(:, :)
    !> Electric field
    real(dp), intent(out)              :: e_field(3)
    !> forces acting on all atoms
    type(force), intent(out)           :: forces


    timeStepMultiplier = int( molecular_dynamics%time_step/timeStepRTTDDFT )
    molecular_dynamics%time_step = timeStepMultiplier*timeStepRTTDDFT
    
    call MD_allocate_global_arrays( nspecies )
    call MD_evaluate_charge_val
    
    call forces%allocate_arrays( natmtot )
    e_field = 0.0_dp
    call force_rttdft( forces, e_field, molecular_dynamics )
    
    allocate( atom_velocities(3, natmtot) )
    call init_atoms_velocities( atom_velocities )

    allocate( atom_positions(3, natmtot) )
    call init_atoms_positions( atom_positions )
    
    if( molecular_dynamics%basis_derivative ) call Update_basis_derivative( atom_velocities, mathcalB, B_time, B_past )
    
    if ( rank == 0 ) then
      call MD_outputs%open_files( natmtot, molecular_dynamics%print_all_force_components  )
      call write_MD_outputs( t_0, atom_positions, atom_velocities, forces, &
                molecular_dynamics%print_all_force_components, MD_outputs )
    end if

  end subroutine

  !> Initialize the velocities of each atom (needed for Ehrenfest MD)
  subroutine init_atoms_velocities(at_velocities)
    !> Velocities of the nuclei at time \( t = 0 \)
    real(dp), intent(inout)         :: at_velocities(:,:)
    integer(i32) :: is, ia, ias

    call assert( size(at_velocities, 1) == 3, 'at_velocities must have size = 3 along dim = 1' )
    call assert( size(at_velocities, 2) == natmtot, 'at_velocities must have size = natmtot along dim = 2' )

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        at_velocities(1:3,ias) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%velocity(1:3)
      end do
    end do
  end subroutine

  !> Initialize the positions of each atom from the global variable `atposc`. 
  !> Needed for Ehrenfest MD
  subroutine init_atoms_positions(at_positions)
    !> Positions of the nuclei at time \( t = 0 \)
    real(dp), intent(inout)         :: at_positions(:,:)
    integer(i32) :: is, ia, ias

    call assert( size(at_positions, 1) == 3, 'at_positions must have size = 3 along dim = 1' )
    call assert( size(at_positions, 2) == natmtot, 'at_positions must have size = natmtot along dim = 2' )
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        at_positions(1:3, ias) = atposc(1:3, ia, is)
      end do
    end do
  end subroutine

  subroutine write_MD_outputs( t, at_positions, at_velocities, forces, print_all_force_components, MD_outputs )
    !> time \(t\)
    real(dp), intent(in) :: t
    !> Positions of the nuclei at time \( t  \)
    real(dp), intent(in)         :: at_positions(:,:)
    !> Velocities of the nuclei at time \( t \)
    real(dp), intent(in)         :: at_velocities(:,:)
    !> forces acting on all atoms
    type(force), intent(in)      :: forces
    !> if `.True.`, print out all contributions to the total force
    logical, intent(in)          :: print_all_force_components
    !> variable with interfaces to MD outputs
    type(MD_out), intent(inout)  :: MD_outputs
    
    if ( print_all_force_components ) then
      call MD_outputs%write_to_files( t, at_positions, at_velocities, forces%total, &
        forces )
    else 
      call MD_outputs%write_to_files( t, at_positions, at_velocities, forces%total )
    end if
  end subroutine

  subroutine deallocate_global_arrays( predictorCorrector, deallocate_ehrenfest_arrays )
    !> when `.True.`, also deallocate arrays used in the predictor corrector loop
    logical, intent(in) :: predictorCorrector
    !> when `.True.`, also deallocate arrays used in Ehrenfest MD
    logical, intent(in) :: deallocate_ehrenfest_arrays

    if ( allocated(apwalm) ) deallocate( apwalm )
    if ( allocated(evecfv_gnd) ) deallocate( evecfv_gnd )
    if ( allocated(evecsv) ) deallocate( evecsv )
    if ( allocated(evecfv_time) ) deallocate( evecfv_time )
    if ( allocated(overlap) ) deallocate( overlap )
    if ( allocated(ham_time) ) deallocate( ham_time )
    if ( allocated(ham_past) ) deallocate( ham_past )
    if ( allocated(pmat) ) deallocate( pmat )
    if ( predictorCorrector ) deallocate( ham_predcorr, evecfv_save )
    if ( nkicks >= 1 ) then
      deallocate( wkick, dirkick, amplkick, t0kick )
    end if
    if ( ntrapcos >= 1 ) then
      deallocate( dirtrapcos, ampltrapcos, omegatrapcos, phasetrapcos )
      deallocate( t0trapcos, trtrapcos, wtrapcos )
    end if
    if ( nsinsq >= 1 ) then
      deallocate( dirsinsq, amplsinsq, omegasinsq )
      deallocate( phasesinsq, t0sinsq, tpulsesinsq )
    end if
    
    if ( deallocate_ehrenfest_arrays ) then
      call MD_deallocate_global_arrays
      if ( allocated( pmatmt ) ) deallocate( pmatmt )
      if ( allocated( mathcalH ) ) deallocate( mathcalH )
      if ( allocated( mathcalB ) ) deallocate( mathcalB )
      if ( allocated( B_past ) ) deallocate( B_past )
      if ( allocated( B_time ) ) deallocate( B_time )
    end if
  end subroutine
end module rttddft_main
