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
!> a RT-TDDFT calculation. 
module rttddft_main
  use asserts, only: assert
  use constants, only: zi
  use MD, only: force, MD_input_keys
  use MD_io, only: MD_out
  use mod_atoms, only: natmtot, natoms, nspecies, atposc, idxas
  use mod_charge_and_moment, only: chgval
  use mod_eigenvalue_occupancy, only: occsv
  use mod_eigensystem, only: nmat
  use mod_kpoint, only: nkpt, wkpt, kpt_latt => vkl
  use mod_lattice, only: omega
  use mod_misc, only: filext
  use mod_mpi_env, only: mpiinfo
  use mod_potential_and_density, only: rhomt, rhoir
  use modinput, only: input, input_type
  use modmpi, only: rank, mpi_env_k, distribute_loop, barrier, terminate_if_false
  use propagators, only: create_propagator, propagator_type => propagator
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_CurrentDensity, only: Current_Density, Current_Density_Field
  use rttddft_Density, only: UpdateDensity
  use rttddft_electric_field, only: Electric_Field, obtain_electric_field
  use rttddft_Energy, only: TotalEnergy, obtain_energy_rttddft
  use rttddft_GlobalVariables
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_init, only: initialize_rttddft
  use rttddft_input, only: rttddft_input_keys
  use rttddft_io, only: open_files_jpa, close_files_jpa, read_jpa, write_jpa, &
    open_file_timing, close_file_timing, write_timing, &
    open_file_nexc, close_file_nexc, write_nexc, &
    open_file_etot, close_file_etot, write_total_energy, &
    open_file_info, close_file_info, write_file_info, write_file_info_header, &
    write_wavefunction, t, RTDDFT_suffix
  use rttddft_MD, only: force_rttdft, move_ions, &
    MD_allocate_global_arrays => allocate_global_arrays, &
    MD_deallocate_global_arrays => deallocate_global_arrays, &
    MD_evaluate_charge_val => evaluate_charge_val
  use rttddft_NumberExcitations, only: Obtain_number_excitations
  use rttddft_Polarization, only: Polarization
  use rttddft_screenshot, only: screenshot
  use rttddft_solve_fields, only: update_a_ind_and_p_vec
  use rttddft_timings, only: Timing_RTTDDFT_and_MD, Timing_RTTDDFT_density, Timing_RTTDDFT_potential, Print_Timings, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential, Vector_Potential_Field
  use rttddft_Wavefunction, only: Update_basis_derivative, normalize_wavefunctions
  use to_char_conversion, only: to_char
  
  implicit none

  private
  
  integer(i32), parameter :: i_spin = 1

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

    ! Basis-expansion coefficients of the groundstate KS-WFs
    ! (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), allocatable :: evecfv_gnd(:, :, :)
    ! Basis-expansion coefficients of the KS-WFs at time \(t\)
    ! (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), allocatable :: evecfv_time(:, :, :) 
    ! Basis-expansion coefficients of the KS-WFs at time \(t\) - auxiliary 
    ! variable used in the predictor-corrector loop
    ! (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), allocatable :: evecfv_save(:, :, :) 
    ! Basis-expansion coefficients of the KS-WFs: second-variational coefficients
    ! (nstfv, nstfv, first_kpt : last_kpt)
    complex(dp), allocatable :: evecsv(:, :, :)
    ! Overlap matrix (of basis functions)
    ! (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), allocatable :: overlap(:, :, :)
    ! Hamiltonian matrix at current time \(t\)
    ! (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), allocatable :: ham_time(:, :, :)
    ! Hamiltonian matrix at previous time \(t - \Delta t \)
    ! (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), allocatable :: ham_past(:, :, :)
    ! Matching coefficients of the (L)APWs
    ! (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), allocatable :: apwalm(:, :, :, :, :)

    ! Momentum matrix elements (projected onto the (L)APW+LO basis elements)
    complex(dp), allocatable  :: pmat(:, :, :, :)
    ! Muffin-tin part of the Momentum matrix
    complex(dp), allocatable  :: pmatmt(:, :, :, :, :)

    integer :: it, first_kpt, last_kpt, n_steps
    integer :: i_print, timeStepMultiplier, l_rad_step
    logical :: predCorrReachedMaxSteps, my_rank_writes_to_output, &
      density_needed, evolve_H0, take_screenshot
    complex(dp), allocatable :: ham_init(:, :, :)
    real(dp), allocatable :: rhoir_init(:), rhomt_init(:, :, :)
    character(len=100)      :: string

    type(Vector_Potential)         :: vec_pot
    type(Vector_Potential_Field)   :: a_ind_save, a_tot_save
    type(Polarization)             :: p_vec, p_vec_save
    type(Current_Density)          :: j_ind, j_ind_save
    ! Spurious paramagnetic current density (obtained for \(t=0\) - this should
    ! ideally be zero for a dense `k-grid` mesh)
    type(Current_Density_Field)    :: j_para_spurious
    type(Electric_Field)           :: e_field

    real(dp), allocatable   :: atom_positions(:, :) ! in cartesian coordinates x, y, z
    real(dp), allocatable   :: atom_velocities(:, :) ! in cartesian coordinates x, y, z
    type(force)             :: forces
    type(MD_input_keys)     :: molecular_dynamics
    type(rttddft_input_keys):: rt
    class(propagator_type), allocatable :: propagator

    ! Current time \( t \) for the time evolution carried out in RT-TDDFT
    real(dp)                :: time

    real(dp), allocatable   :: n_exc(:), n_gs(:)
    real(dp)                :: timei, timef, timeiter, dt
    real(dp)                :: tol, eps_occ
    real(dp), parameter     :: tol_default = 1e-10_dp
    type(MD_out)            :: MD_outputs

    ! Variables to store data and print
    real(dp), allocatable    :: time_store(:)
    type(Vector_Potential_Field), allocatable :: a_ind_store(:), a_tot_store(:)
    type(Current_Density_Field), allocatable  :: j_ind_store(:)
    type(Polarization), allocatable :: p_vec_store(:)
    real(dp),allocatable    :: atposcstore(:,:,:), velstore(:,:,:)
    type(force),allocatable :: forces_store(:)
    logical,allocatable     :: print_forces(:)

    type(TotalEnergy), allocatable  :: etotstore(:)
    type(Timing_RTTDDFT_and_MD)     :: timing
    type(Timing_RTTDDFT_and_MD), allocatable :: timing_store(:)

    call timesec( timei )

    ! Sanity check
    call sanity_checks( input )

    ! Interface with input parameters
    tol = tol_default
    if( associated(input%groundstate%solver) ) tol = input%groundstate%solver%evaltol
    call rt%parse_input( input%xs%realTimeTDDFT, tol, vec_pot )
    call molecular_dynamics%parse_input()
    ! we only perform MD in RT-TDDFT if the type is Ehrenfest
    if( molecular_dynamics%on ) molecular_dynamics%on = ( trim(molecular_dynamics%MD_type) == 'Ehrenfest' )
    call create_propagator( propagator, rt%propagator_input, .not. molecular_dynamics%on )
    
    ! Output general info to RTTDDFT_INFO.OUT
    my_rank_writes_to_output = ( rank == 0 ) 
    if( my_rank_writes_to_output ) then
      call open_rttddft_outputs( rt )
      call write_file_info_header
    end if

    ! Initialization
    time = 0._dp
    dt = rt%propagator_input%dt()
    n_steps = int( rt%t_end / dt )
    l_rad_step = input%groundstate%lradstep
    eps_occ = input%groundstate%epsocc
    call initialize_rttddft( rt%pmat, rt%predictor_corrector%on, propagator%extrapolation_needed(), &
      vec_pot, molecular_dynamics, evecfv_gnd, evecfv_time, evecfv_save, evecsv, &
      overlap, ham_time, ham_past, apwalm, pmat, pmatmt )
    if( molecular_dynamics%on ) call init_MD( time, vec_pot%a_tot, dt, &
      evecfv_time, overlap, ham_time, timeStepMultiplier, molecular_dynamics, &
      MD_outputs, atom_positions, atom_velocities, e_field, forces )
    if ( rt%subtract_J0 ) then
      call j_ind%evaluate_paramagnetic( evecfv_gnd, pmat, occsv(:, first_kpt:last_kpt), wkpt(first_kpt:last_kpt), mpi_env_k )
      j_para_spurious = j_ind%paramagnetic
    end if
    
    ! Allocate variables to be stored and printed only after rt_input%n_print steps
    allocate( time_store(rt%n_print), a_ind_store(rt%n_print), a_tot_store(rt%n_print))
    allocate( j_ind_store(rt%n_print), p_vec_store(rt%n_print) )
    if( rt%printTimings%general() ) then
      allocate( timing_store(rt%n_print) )
    end if
    if( rt%calculate_total_energy ) allocate(etotstore(rt%n_print))
    if( rt%calculate_n_exc ) allocate( n_exc(rt%n_print), n_gs(rt%n_print) )
    if( molecular_dynamics%on ) then
      allocate( print_forces(rt%n_print), atposcstore(3,natmtot,rt%n_print), velstore(3,natmtot,rt%n_print))
      if ( molecular_dynamics%print_all_force_components ) then
        allocate( forces_store(rt%n_print) )
        do i_print = 1, rt%n_print
          call forces_store(i_print)%allocate_arrays( natmtot )
        end do
      end if
    end if

    if( my_rank_writes_to_output ) call write_fields( [time], [vec_pot%a_ind], [vec_pot%a_tot], & 
      [p_vec], [j_ind%total()] )

    ! Initialize integers that contain the first and last k-point
    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    ! Total energy
    if ( rt%calculate_total_energy ) then
      call potcoul
      call potxc
      call obtain_energy_rttddft( first_kpt, ham_time, evecfv_gnd, mpi_env_k, etotstore(1) )
      if( my_rank_writes_to_output ) call write_total_energy( .True., [time], [etotstore(1)] )
    end if

    ! Number of excitations
    if (rt%calculate_n_exc) then
      call Obtain_number_excitations( evecfv_gnd, evecfv_time, overlap, eps_occ, &
        & occsv(:, first_kpt:last_kpt), wkpt(first_kpt:last_kpt), mpi_env_k, n_exc(1), n_gs(1) )
      if( my_rank_writes_to_output ) call write_nexc( .True., [time], [n_exc(1)], [n_gs(1)] )
    end if

    if ( rt%screenshots%on ) then
      if ( rt%screenshots%density%on ) then
        call UpdateDensity( first_kpt, evecfv_time(:, :, first_kpt : last_kpt), &
          evecsv, it, rt%normalize_WF, l_rad_step, rt%printTimings, timing%t_RTTDDFT%dens )
        rhomt_init = rhomt
        rhoir_init = rhoir
      end if
      call screenshot( 0, rt%screenshots, overlap, evecfv_gnd, evecfv_time, ham_time, nmat(1, first_kpt:last_kpt), occsv(:, first_kpt:last_kpt), rhomt, rhoir, mpi_env=mpi_env_k )
    end if ! rt%screenshots%on

    if( rt%printTimings%general() ) then
      call timesec( timef )
      if( my_rank_writes_to_output ) call write_timing( timef-timei ) !write time for initialization
    end if


    ! whether explicitly field-independent Hamiltonian should be evolved in time
    evolve_H0 = .false.
    if ( molecular_dynamics%on .or. ( .not. rt%eeInteraction%ipa ) ) evolve_H0 = .true.
    if ( .not. evolve_H0  ) allocate( ham_init, source = ham_time )

    i_print = 1
    timeiter = timef
    ! This is the most important loop (performed for each time step \(\Delta t\)
    do it = 1, n_steps

      call timing%reset()
      ! Variable to store the timing of each iteration
      timei = timeiter

      ! The "real time" t of our evolution
      time = time + dt

      ! Shall the screenshot be taken on the current step
      take_screenshot = .false.
      if ( rt%screenshots%on ) take_screenshot = ( mod( it, rt%screenshots%n_steps ) == 0 )
        
      ! We may need to update charge density on some steps
      density_needed = ( .not. rt%eeInteraction%ipa )
      if ( take_screenshot ) density_needed = density_needed .or. rt%screenshots%density%on

      ! WAVEFUNCTION
      if ( rt%predictor_corrector%on ) evecfv_save = evecfv_time
      if ( molecular_dynamics%on .and. molecular_dynamics%basis_derivative ) then
        call Update_basis_derivative( atom_velocities, mathcalB, B_time, B_past )
        ham_time = ham_time - zi*B_time
      end if
      call propagator%evolve( list_of_H_minus_dt=ham_past, list_of_H_0=ham_time, list_of_S=overlap, psi=evecfv_time, dims=nmat(i_spin, first_kpt:last_kpt) )
      if ( rt%normalize_WF ) call normalize_wavefunctions( overlap, evecfv_time )
      if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%wavefunction )

      ! Update the paramagnetic component of the induced current density
      j_ind_save = j_ind
      call j_ind%evaluate_paramagnetic( evecfv_time, pmat, occsv(:, first_kpt:last_kpt), &
        wkpt(first_kpt:last_kpt), mpi_env_k )
      if ( rt%subtract_J0 ) call j_ind%paramagnetic%add_vector( -j_para_spurious%components )
      if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%current_density )

      ! DENSITY
      if ( density_needed ) call UpdateDensity( first_kpt, evecfv_time(:, :, first_kpt : last_kpt), &
        evecsv, it, rt%normalize_WF, l_rad_step, rt%printTimings, timing%t_RTTDDFT%dens )

      ! KS-POTENTIAL
      if ( .not. rt%eeInteraction%ipa ) call uppot( rt%printTimings, timing%t_RTTDDFT%pot )

      ! VECTOR POTENTIAL
      if( rt%printTimings%general() ) call timesec( timei )
      ! Check if we need to save aind, pvec, atot and aext
      if( vec_pot%is_external_field_given() .and. rt%predictor_corrector%on .and. ( .not. vec_pot%is_solver_euler() ) ) then
        a_ind_save = vec_pot%a_ind
        p_vec_save = p_vec
      end if
      a_tot_save = vec_pot%a_tot
      call update_a_ind_and_p_vec( time, dt, j_ind_save, j_ind%paramagnetic, vec_pot, p_vec )
      call vec_pot%evaluate_a_tot( time )
      if( molecular_dynamics%on ) then
        if( vec_pot%is_total_field_given() ) then
          call e_field%obtain_electric_field( 2*dt, Vector_Potential_Field(vec_pot%applied_vector_potential( time+dt )), a_tot_save )
        else
          call e_field%obtain_electric_field( dt, vec_pot%a_tot, a_tot_save )
        end if
      end if
      if( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%vector_potential )

      ! INDUCED CURRENT
      ! Update the diamagnetic component of the induced current density
      call j_ind%evaluate_diamagnetic( chgval/Omega, vec_pot%a_tot )
      ! Update the total induced current
      if ( rt%predictor_corrector%on .and. ( .not. vec_pot%is_solver_euler() ) ) j_ind_save = j_ind

      ! HAMILTONIAN
      if( propagator%extrapolation_needed() ) ham_past = ham_time
      call UpdateHam( first_kpt, vec_pot%a_tot, calculateOverlap=.False., &
        calculateH0=evolve_H0, forcePmatHermitian=rt%pmat%force_pmat_hermitian, &
        overlap=overlap, ham_time=ham_time, apwalm=apwalm, pmat=pmat, &
        printTimings=rt%printTimings, t_ham=timing%t_RTTDDFT%ham, t_MD=timing%t_Ehrenfest, &
        update_mathcalH=.False., update_mathcalB=.False., update_pmat=.False., ham_init=ham_init )

      if ( rt%predictor_corrector%on ) then
        if ( rt%printTimings%general() ) call timesec( timei )
        call loopPredictorCorrector( it, time, rt, l_rad_step, first_kpt, &
          evecfv_time, evecfv_save, evecsv, overlap, ham_time, ham_past, apwalm, pmat, &
          a_ind_save, a_tot_save, p_vec_save, j_ind_save, j_para_spurious, &
          propagator, vec_pot, p_vec, j_ind, mpi_env_k, predCorrReachedMaxSteps )
        if ( predCorrReachedMaxSteps .and. my_rank_writes_to_output ) call warning( 'Problems with convergence (PredCorr), time: ' //  to_char(time) )
        if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%pred_corr )
      end if !predictor-corrector

      ! Obtain the total energy, if requested
      if( rt%calculate_total_energy ) then
        if ( rt%printTimings%detailed() ) call timesec( timei )
        call obtain_energy_rttddft( first_kpt, ham_time, evecfv_time, mpi_env_k, etotstore(i_print) )
        if ( rt%printTimings%detailed() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%energy )
      end if

      ! Obtain the number of excited electrons, if requested
      if( rt%calculate_n_exc ) then
        if ( rt%printTimings%detailed() ) call timesec( timei )
        call Obtain_number_excitations( evecfv_gnd, evecfv_time, overlap, eps_occ, &
          & occsv(:, first_kpt:last_kpt), wkpt(first_kpt:last_kpt), mpi_env_k, n_exc(i_print), n_gs(i_print))
        if( rt%printTimings%detailed() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%n_exc )
      end if

      if ( molecular_dynamics%on ) then
        print_forces(i_print) = .False.
        if ( mod( it, timeStepMultiplier ) == 0 ) then
          if ( rt%printTimings%general() ) call timesec( timei )
          call forces%save_total_force()
          call force_rttdft( forces, vec_pot%a_tot, e_field, molecular_dynamics, evecfv_time, overlap, ham_time, rt%printTimings, timing%t_Ehrenfest )
          call move_ions( first_kpt, forces%total, forces%total_save, molecular_dynamics%time_step, &
            atom_velocities, apwalm, rt%printTimings, timing%t_Ehrenfest )
          print_forces(i_print) = .True.
          call get_atoms_positions( atposcstore(:, :, i_print) )
          velstore(:, :, i_print) = atom_velocities
          forces_store(i_print) = forces
          ! Update Hamiltonian with the new basis
          if( molecular_dynamics%update_overlap .or. allocated(mathcalH) .or. &
            & allocated(mathcalB) .or. molecular_dynamics%update_pmat ) then
            if( propagator%extrapolation_needed() ) ham_past = ham_time  
            call UpdateHam( first_kpt, vec_pot%a_tot, &
              & forcePmatHermitian=rt%pmat%force_pmat_hermitian, &
              & calculateOverlap=molecular_dynamics%update_overlap, calculateH0=evolve_H0, &
              & overlap=overlap, ham_time=ham_time, apwalm=apwalm, pmat=pmat, pmatmt=pmatmt, &
              & printTimings=rt%printTimings, t_ham=timing%t_RTTDDFT%ham, t_MD=timing%t_Ehrenfest, &
              & update_mathcalH=allocated(mathcalH), &
              & update_mathcalB=allocated(mathcalB), &
              & update_pmat=molecular_dynamics%update_pmat, ham_init=ham_init )
          end if
          if( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_Ehrenfest%t_MD_step )
        end if ! if ( mod( it, timeStepMultiplier ) == 0 )
      end if ! if ( molecular_dynamics%on ) then

      ! Check if a screenshot has been requested
      if ( take_screenshot ) then
        if( rt%printTimings%general() ) call timesec( timei )
        call screenshot( it, rt%screenshots, overlap, evecfv_gnd, evecfv_time, ham_time, nmat(1, first_kpt:last_kpt), &
          occsv(:, first_kpt:last_kpt), rhomt, rhoir, rhomt_init, rhoir_init, mpi_env_k )
        if( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%screenshot )
      end if

      ! Store relevant information from this iteration
      time_store(i_print) = time
      a_ind_store(i_print) = vec_pot%a_ind
      a_tot_store(i_print) = vec_pot%a_tot
      p_vec_store(i_print) = p_vec
      j_ind_store(i_print) = j_ind%total()
      if( rt%printTimings%general() ) timing_store(i_print) = timing

      ! Print relevant information, every 'rt%n_print' steps
      if ( i_print == rt%n_print ) then
        if( my_rank_writes_to_output ) then
          call write_fields( time_store, a_ind_store, a_tot_store, p_vec_store, j_ind_store )
          if ( rt%calculate_total_energy ) call write_total_energy( .False., time_store, etotstore )
          if ( rt%calculate_n_exc ) call write_nexc( .False., time_store, n_exc, n_gs )

          if( molecular_dynamics%on ) then
            do i_print = 1, rt%n_print
              if( print_forces(i_print) ) call write_MD_outputs( time_store(i_print), &
                atposcstore(:, :, i_print), velstore(:,:,i_print), forces_store(i_print), &
                molecular_dynamics%print_all_force_components, MD_outputs )
            end do
          end if ! if( molecular_dynamics%on )
        end if

        if( rt%printTimings%general() ) then
          call timesec_RTTDDFT( timeiter, timing_store(rt%n_print)%t_iteration )
          if( my_rank_writes_to_output ) call write_timing( it, timing_store, &
          molecular_dynamics%on )
        end if
        i_print = 0
      else ! if ( iprint == rt%n_print ) 
        if( rt%printTimings%general() ) call timesec_RTTDDFT( timeiter, timing_store(i_print)%t_iteration )
      end if ! if ( iprint == rt%n_print ) 
      i_print = i_print + 1
      ! Make all the processes wait here: the master alone has been writing the files above
      call barrier( mpi_env_k )
    end do ! do it = 1, nsteps

    if ( my_rank_writes_to_output ) then
      call write_file_info( 'Real-time TDDFT calculation finished' )
      call close_rttddft_outputs( rt )
      if ( molecular_dynamics%on ) call MD_outputs%close_files()
    end if

    ! write wavefunction, and potential and density with _RTTDDFT.OUT as suffix
    call write_wavefunction( t, first_kpt, kpt_latt(:, first_kpt:last_kpt), evecfv_time, mpi_env_k )
    string = filext
    filext = RTDDFT_suffix // trim( filext )
    if ( my_rank_writes_to_output ) call writestate
    filext = string

    call deallocate_global_arrays( molecular_dynamics%on )
    
  end subroutine coordinate_rttddft_calculation


  !> Open output files
  subroutine open_rttddft_outputs( rt_input )
    !> Type that encapsulates the input keywords
    type(rttddft_input_keys), intent(in) :: rt_input

    call open_files_jpa( new=.True. )
    call open_file_info
    if( rt_input%calculate_total_energy ) call open_file_etot( new=.True. )
    if( rt_input%calculate_n_exc ) call open_file_nexc( new=.True. )
    if( rt_input%printTimings%general() ) call open_file_timing( new=.True.)

  end subroutine

  !> Close output files
  subroutine close_rttddft_outputs( rt_input )
    !> Type that encapsulates the input keywords
    type(rttddft_input_keys), intent(in) :: rt_input

    call close_files_jpa
    call close_file_info
    if( rt_input%calculate_total_energy ) call close_file_etot
    if( rt_input%calculate_n_exc ) call close_file_nexc
    if( rt_input%printTimings%general() ) call close_file_timing

  end subroutine

  !> (private) Read time and fields stored in the corresponding files
  subroutine read_time_and_fields( t, p_vec, a_t, a_ind_t_minus_dt, a_tot_t_minus_dt )
    !> Time \(t\)
    real(dp), intent(out) :: t
    !> Polarization vector
    type(Polarization), intent(out) :: p_vec
    !> Vector potential at time \(t\)
    type(Vector_Potential), intent(inout) :: a_t
    !> \(\mathbf{A}_{ind}) at time \(t-\Delta t\)
    type(Vector_Potential_Field), intent(out) :: a_ind_t_minus_dt
    !> \(\mathbf{A}_{tot}) at time \(t-\Delta t\)
    type(Vector_Potential_Field), intent(out) :: a_tot_t_minus_dt

    call read_jpa( t, p_vec )
    call read_jpa( t, a_t%a_ind, a_ind_t_minus_dt, a_t%a_tot, a_tot_t_minus_dt)
  end subroutine

  !> Wrapper to call [[write_jpa]]
  subroutine write_fields( time_array, a_ind_array, a_tot_array, p_vec_array, j_ind_array )
    !> Array with times
    real(dp), intent(in) :: time_array(:)
    !> Array with the induced vector fields
    type(Vector_Potential_Field), intent(in) :: a_ind_array(:)
    !> Array with the total vector fields
    type(Vector_Potential_Field), intent(in) :: a_tot_array(:)
    !> Array with the polarization fields
    type(Polarization), intent(in) :: p_vec_array(:)
    !> Array with the current density field
    type(Current_Density_Field), intent(in) :: j_ind_array(:)

    call write_jpa( time_array, a_ind_array, a_tot_array )
    call write_jpa( time_array, p_vec_array )
    call write_jpa( time_array, j_ind_array )
  end subroutine

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

    if( associated(inp%xs%realTimeTDDFT%predictorCorrector) ) then
      ! Consistency check: MD and predictor corrector?
      call terminate_if_false( .not. associated(inp%MD), &
        & 'It is currently not possible to use the predictor corrector method together with molecular dynamics')
      ! Consistency check: predictor corrector method cannot be used with propagators SE and EH
      call terminate_if_false( trim(inp%xs%realTimeTDDFT%propagator)/='SE' .and. trim(inp%xs%realTimeTDDFT%propagator)/='EH', &
        & 'EH and SE methods are not compatible with predictor-corrector' )
      ! Consistency check: predictor corrector method should not be used with frozen ee interaction
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%eeInteraction ) /= "IPA", &
        & 'Predictor corrector method should not be used together with IP approximation')
    end if

    if ( inp%xs%realTimeTDDFT%calculateTotalEnergy ) then
      ! Consistency check: real-time total energy is ill-defined with frozen ee interaction
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%eeInteraction ) /= "IPA", &
        & 'Real-time total energy should not be eveluated with IP approximation')
    end if

  end subroutine

  !> Loop used in the predictor-corrector method
  subroutine loopPredictorCorrector( it, time, rt, l_rad_step, first_kpt, &
    evecfv_time, evecfv_save, evecsv, overlap, ham_time, ham_past, apwalm, pmat, &
    a_ind_t_minus_dt, a_tot_t_minus_dt, p_vec_t_minus_dt, j_t_minus_dt, j_para_spurious,&
    propagator, a_t, p_vec, j_t, mpi_env, maxStepsReached )
    !> current iteration number in the RT-TDDFT loop
    integer(i32), intent(in) :: it
    !> time \( t \)
    real(dp), intent(in) :: time
    !> Type that encapsulates the parameters defined in the input file
    type(rttddft_input_keys), intent(in) :: rt
    !> radial step length
    integer(i32), intent(in) :: l_rad_step
    !> index of the first `k-point` to be considered in the sum
    integer(i32),intent(in) :: first_kpt
    !> Basis-expansion coefficients of the KS-WFs at time \(t\)
    !> (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), intent(out) :: evecfv_time(:, :, first_kpt:)
    !> Basis-expansion coefficients of the KS-WFs at time \(t\) - auxiliary 
    !> variable used in the predictor-corrector loop
    !> (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), intent(in) :: evecfv_save(:, :, first_kpt:)
    !> Basis-expansion coefficients of the KS-WFs: second-variational coefficients
    !> (nstfv, nstfv, first_kpt : last_kpt)
    complex(dp), intent(in) :: evecsv(:, :, first_kpt:)
    !> Overlap matrix (of basis functions)
    !> (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), intent(inout) :: overlap(:, :, first_kpt:)
    !> Hamiltonian matrix at current time \(t\)
    !> (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), intent(inout) :: ham_time(:, :, first_kpt:)
    !> Hamiltonian matrix at previous time \(t - \Delta t \)
    !> (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), intent(inout) :: ham_past(:, :, first_kpt:)
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), intent(in) :: apwalm(:, :, :, :, first_kpt:)
    !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
    complex(dp), intent(inout) :: pmat(:, :, :, first_kpt:)
    !> `aind` at time \( t-\Delta t\) 
    class(Vector_Potential_Field), intent(in) :: a_ind_t_minus_dt
    !> `atot` at time \( t-\Delta t\) 
    class(Vector_Potential_Field), intent(in) :: a_tot_t_minus_dt
    !> Polarization at time \( t-\Delta t\) 
    type(Polarization) :: p_vec_t_minus_dt
    !> Current density at time \( t-\Delta t\) 
    class(Current_Density), intent(in) :: j_t_minus_dt
    !> Spurious paramagnetic current density (obtained at \(t=0\))
    class(Current_Density_Field), intent(in) :: j_para_spurious
    !> Propagator
    class(propagator_type), intent(in) :: propagator
    !> Structure with the vector potential
    type(Vector_Potential), intent(inout) :: a_t
    !> Polarization at time \( t \)
    type(Polarization), intent(inout) :: p_vec
    !> Current density at time \( t) 
    type(Current_Density), intent(inout) :: j_t
    !> MPI environment
    type(mpiinfo), intent(in)  :: mpi_env
    !> When `.True.`, it informs that the maximum steps have been reached
    logical, intent(out) :: maxStepsReached

    integer(i32) :: i, last_kpt, nham
    real(dp)     :: err, dt
    complex(dp), allocatable :: ham_predcorr(:, :, :)

    dt = rt%propagator_input%dt()
    last_kpt = ubound( evecfv_time, 3 )
    nham = size( ham_time, 1 )
    allocate( ham_predcorr(nham, nham, first_kpt:last_kpt) )

    do i = 1, rt%predictor_corrector%max_steps
      ! WAVEFUNCTION
      evecfv_time = evecfv_save
      call propagator%evolve( list_of_H_dt=ham_time, list_of_H_0=ham_past, list_of_S=overlap, psi=evecfv_time, dims=nmat(i_spin, first_kpt:last_kpt) )
      if ( rt%normalize_WF ) call normalize_wavefunctions( overlap, evecfv_time )

      ! Update the paramagnetic component of the induced current density
      j_t = j_t_minus_dt
      call j_t%evaluate_paramagnetic( evecfv_time, pmat, occsv(:, first_kpt:last_kpt), &
        wkpt(first_kpt:last_kpt), mpi_env )
      if ( rt%subtract_J0 ) call j_t%paramagnetic%add_vector( -j_para_spurious%components )

      ! DENSITY
      call UpdateDensity( first_kpt, evecfv_time, evecsv, it, rt%normalize_WF, l_rad_step )
      ! KS-POTENTIAL
      call uppot()

      ! VECTOR POTENTIAL
      ! Update the induced part of the vector potential
      if( a_t%is_external_field_given() ) then
        call a_t%set_a_tot_a_ind( a_tot_t_minus_dt, a_ind_t_minus_dt )
        p_vec = p_vec_t_minus_dt
        call update_a_ind_and_p_vec( time, dt, j_t_minus_dt, j_t%paramagnetic, a_t, p_vec )
        call a_t%evaluate_a_tot( time )
      end if

      ! INDUCED CURRENT
      ! Update the paramagnetic component of the induced current density
      call j_t%evaluate_diamagnetic( chgval/Omega, a_t%a_tot )

      ! HAMILTONIAN
      ham_predcorr = ham_time
      ham_time = ham_past
      call UpdateHam( first_kpt, a_t%a_tot, calculateOverlap=.False., &
        calculateH0=.true., forcePmatHermitian=rt%pmat%force_pmat_hermitian, &
        overlap=overlap, ham_time=ham_time, apwalm=apwalm, &
        pmat=pmat )

      ! Check the difference between the two hamiltonians
      err = maxval( abs(ham_predcorr - ham_time) )
      if ( err <= rt%predictor_corrector%tol ) exit

    end do
    maxStepsReached = (i>rt%predictor_corrector%max_steps)
  end subroutine 

  !> Subroutine to initialize all MD related variables
  subroutine init_MD( t_0, a_tot, timeStepRTTDDFT, evecfv_time, overlap, ham_time, timeStepMultiplier, &
    molecular_dynamics, MD_outputs, atom_positions, atom_velocities, e_field, forces )
    !> Initial time \( t_0 \)
    real(dp), intent(in)               :: t_0
    !> Vector potential (total)
    class(Vector_Potential_Field), intent(in) :: a_tot
    !> Time step used in the real-time TDDFT calculation
    real(dp), intent(in)               :: timeStepRTTDDFT
    !> Basis-expansion coefficients of the KS-WFs at time \(t\)
    complex(dp), intent(in) :: evecfv_time(:, :, :)
        !> Overlap matrix (of basis functions)
    complex(dp), intent(in) :: overlap(:, :, :)
    !> Hamiltonian matrix at current time \(t\)
    complex(dp), intent(in) :: ham_time(:, :, :)
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
    type(Electric_Field), intent(in)   :: e_field
    !> forces acting on all atoms
    type(force), intent(out)           :: forces


    timeStepMultiplier = int( molecular_dynamics%time_step/timeStepRTTDDFT )
    molecular_dynamics%time_step = timeStepMultiplier*timeStepRTTDDFT
    
    call MD_allocate_global_arrays( nspecies )
    call MD_evaluate_charge_val
    
    call forces%allocate_arrays( natmtot )
    call force_rttdft( forces, a_tot, e_field, molecular_dynamics, evecfv_time, &
      overlap, ham_time )
    
    allocate( atom_velocities(3, natmtot) )
    call init_atoms_velocities( atom_velocities )

    allocate( atom_positions(3, natmtot) )
    call get_atoms_positions( atom_positions )
    
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

  !> Obtain the positions of each atom from the global variable `atposc`. 
  !> Needed for Ehrenfest MD
  subroutine get_atoms_positions(at_positions)
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

  subroutine deallocate_global_arrays( deallocate_ehrenfest_arrays )
    !> when `.True.`, also deallocate arrays used in Ehrenfest MD
    logical, intent(in) :: deallocate_ehrenfest_arrays
    
    if ( deallocate_ehrenfest_arrays ) then
      call MD_deallocate_global_arrays
      if ( allocated( mathcalH ) ) deallocate( mathcalH )
      if ( allocated( mathcalB ) ) deallocate( mathcalB )
      if ( allocated( B_past ) ) deallocate( B_past )
      if ( allocated( B_time ) ) deallocate( B_time )
    end if
  end subroutine
end module rttddft_main
