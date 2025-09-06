! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module is the kernel of a RT-TDDFT calculation.
!> It contains the subroutine `coordinate_rttddft_calculation`, which manages
!> a RT-TDDFT calculation. 
module rttddft_main
  use asserts, only: assert
  use constants, only: zi, real_zero, zzero
  use matrix_elements, only: me_finit
  use MD, only: force, MD_input_keys, trajectory
  use MD_io, only: MD_out
  use mod_atoms, only: natmtot, natoms, nspecies, atposc, idxas
  use mod_charge_and_moment, only: chgval
  use mod_kpointset, only: Gk_set, G_set, k_set
  use mod_lattice, only: omega, avec
  use mod_mpi_env, only: mpiinfo
  use mod_potential_and_density, only: rhomt, rhoir
  use modinput, only: input, input_type
  use modmpi, only: mpiglobal, mpi_env_k, distribute_loop, barrier, terminate_if_false
  use physical_constants, only: c
  use precision, only: dp, i32, sp
  use propagators, only: propagator_type => propagator
  use rttddft_berry, only: get_td_overlap_det_and_berry_coupling_term
  use rttddft_CurrentDensity, only: Current_Density, Current_Density_Field
  use rttddft_Density, only: update_density, ground_state
  use rttddft_file_names, only: RTTDDFT_suffix
  use rttddft_electric_field, only: Electric_Field, obtain_electric_field
  use rttddft_Energy, only: TotalEnergy, obtain_energy_rttddft
  use rttddft_GlobalMDVariables
  use rttddft_Hamiltonian, only: add_external_coupling_berry_phase, add_external_coupling_velocity_gauge, &
    update_hamiltonian_without_pa_term_ks, update_hamiltonian_without_pa_term_lapw
  use rttddft_init, only: initialize_me, initialize_rttddft
  use rttddft_input, only: rttddft_input_keys
  use rttddft_io, only: open_files_vector_fields, close_files_vector_fields, read_vector_field, write_vector_field, &
    open_file_timing, close_file_timing, write_timing, &
    open_file_nexc, close_file_nexc, write_nexc, &
    open_file_etot, close_file_etot, write_total_energy, &
    open_file_info, close_file_info, write_file_info, write_file_info_header, &
    write_wavefunction, t, t_minus_dt, copy_files, write_state_Ehrenfest_MD, read_state_Ehrenfest_MD, &
    write_phases
  use rttddft_MD, only: force_rttdft, move_ions, update_basis_derivative, &
    MD_allocate_global_arrays => allocate_global_arrays, &
    MD_deallocate_global_arrays => deallocate_global_arrays, &
    MD_evaluate_charge_val => evaluate_charge_val
  use rttddft_Overlap, only: overlap_set
  use rttddft_Polarization, only: Polarization
  use rttddft_pmat, only: obtain_pmat_LAPWLOBasis
  use rttddft_potential, only: update_potential
  use rttddft_sanity_checks, only: check_rttddft_input, check_rttddft_setup
  use rttddft_screenshot, only: screenshot
  use rttddft_solve_fields, only: update_a_ind_and_p_vec
  use rttddft_timings, only: Timing_RTTDDFT_and_MD, Timing_RTTDDFT_density, &
    Timing_RTTDDFT_potential, Print_Timings, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential, Vector_Potential_Field
  use rttddft_Wavefunction, only: wavefunction_set
  use to_char_conversion, only: to_char

  implicit none

  private
  
  integer(i32), parameter :: i_spin = 1

  public :: coordinate_rttddft_calculation

contains

  !> This subroutine manages a RT-TDDFT calculation.  
  !> 1. Run a single-shot groundstate calculation using the already converged
  !> density and potential.
  !> 2. Obtain the KS wavefunctions, the density and the hamiltonian at 
  !> \( t=0 \), from the previous step.
  !> 3. If restarting from a previous calculation (`do="fromfile"`), 
  !> then, using the results from the previous calculation: 
  !>    1. set the current time \(t\), and the vector fields using the results; 
  !>    2. read the wavefunctions at time \(t\), and, if needed, at time \(t-\Delta t\);
  !>    3. using these wavefunctions, determine the density, the potential, and the hamiltonian.
  !> 4. Evolve the wavefunctions, the density and the hamiltonian using the desired time step.
  subroutine coordinate_rttddft_calculation()
    ! Basis-expansion coefficients of the KS-WFs
    class(wavefunction_set), allocatable :: psi
    ! Overlap and Hamiltonian matrices (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), allocatable :: ham_time(:, :, :), berry_coupling_term(:, :, :), &
      ham_past(:, :, :), ham_init(:, :, :), effective_potential_init(:, :, :)
    type(overlap_set) :: overlap
    ! Matching coefficients of the (L)APWs: (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), allocatable :: apwalm(:, :, :, :, :)
    ! Momentum matrix elements (ham_dimension, ham_dimension, 3, first_kpt : last_kpt)
    ! and (nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt)
    complex(dp), allocatable :: pmat(:, :, :, :), pmatmt(:, :, :, :, :)
    ! Electron density (lmmaxvr, nrmtmax, natmtot) and (ngrtot)
    real(dp), allocatable :: rhomt_frozen(:, :, :), rhomt_init(:, :, :), &
      rhoir_frozen(:), rhoir_init(:)
    ! Initial occupations and energies array (nstates, first_kpt : last_kpt)
    real(dp), allocatable :: occupations(:, :), initial_ks_energies(:, :)
    ! k-dependent Hamiltonian's dimensions array (first_kpt : last_kpt)
    integer(i32), allocatable :: k_dependent_dims(:)
    ! KS-LAPW+lo transition matrix (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), allocatable :: ks_lapwlo_transition_matrix(:, :, :)
    ! Planewave matrix elements between neighbouring k points
    complex(dp), allocatable :: pws_for_berry_phase(:, :, :, :, :)
    ! Indices of the neighbouring k points
    integer(i32), allocatable :: k_ptrs(:, :, :)

    integer(i32) :: it, first_kpt, last_kpt, first_step, last_step, i_print, &
      time_step_multiplier, lmax_potential
    logical :: pred_corr_reached_max_steps, my_rank_writes_to_output, &
      density_needed, evolve_H0, take_screenshot
    type(Vector_Potential) :: vec_pot
    type(Vector_Potential_Field) :: a_ind_save, a_tot_save
    type(Polarization) :: p_vec, p_vec_save, p_vec_init, p_vec_prev
    type(Electric_Field) :: e_field, e_vec_save
    type(Current_Density) :: j_ind, j_ind_save
    ! Spurious paramagnetic current density (obtained at \(t=0\) - it should ideally be zero for a dense `k-grid` mesh)
    type(Current_Density_Field) :: j_para_spurious
    
    type(trajectory) :: nuclei_motion
    type(force) :: forces
    type(MD_input_keys) :: molecular_dynamics
    type(rttddft_input_keys) :: rt
    class(propagator_type), allocatable :: propagator

    type(k_set) :: kset
    type(Gk_set) :: Gkset
    ! Note: Gset is not explicitly used, but me_basis generated for ME
    ! evaluation points to it implicitly
    type(G_set) :: Gset
    complex(dp), allocatable :: td_overlap_det(:, :)

    real(dp) :: time, timei, timef, time_aux, timeiter, dt, tol, eps_occ, energy_gap
    real(dp), allocatable :: n_exc(:), n_gs(:), prev_phases(:, :)
    real(dp), parameter :: tol_default = 1e-10_dp
    type(MD_out) :: MD_outputs
    ! Variables to store data and print
    real(dp), allocatable :: time_store(:)
    type(Vector_Potential_Field), allocatable :: a_ind_store(:), a_tot_store(:)
    type(Electric_Field), allocatable :: e_vec_store(:)
    type(Current_Density_Field), allocatable :: j_ind_store(:)
    type(Polarization), allocatable :: p_vec_store(:)
    type(force), allocatable :: forces_store(:)
    type(trajectory), allocatable :: nuclei_motion_store(:)
    logical, allocatable :: print_forces(:)
    type(TotalEnergy), allocatable :: etotstore(:)
    type(Timing_RTTDDFT_and_MD) :: timing
    type(Timing_RTTDDFT_and_MD), allocatable :: timing_store(:)

    call timesec( timei )

    ! Input sanity check
    call check_rttddft_input( input )

    ! Interface with input parameters
    tol = tol_default
    if( associated( input%groundstate%solver ) ) tol = input%groundstate%solver%evaltol
    call rt%parse_input( input, tol, vec_pot )
    call molecular_dynamics%parse_input()
    lmax_potential = input%groundstate%lmaxvr
    ! we only perform MD in RT-TDDFT if the type is Ehrenfest
    if( molecular_dynamics%on ) molecular_dynamics%on = ( trim(molecular_dynamics%MD_type) == 'Ehrenfest' )
    
    my_rank_writes_to_output = mpiglobal%is_root 
    if( rt%restart_previous_calculation() ) then
      if( rt%restart_extension /= "" ) then
        ! Copy files: sources are files ending with `rt%restart_extension`, dest. are to the default file names
        if( my_rank_writes_to_output ) then 
          call copy_files( rt%restart_extension, rt%calculate_n_exc, rt%calculate_total_energy )
          if( molecular_dynamics%on ) call MD_outputs%copy_files( rt%restart_extension, molecular_dynamics%print_all_force_components )
        end if
        ! Before reading, ensure that copying has been finished
        call barrier()
      end if
      call read_time_and_fields( time, p_vec, vec_pot, a_ind_save, a_tot_save, e_field, e_vec_save )
      if( molecular_dynamics%on ) then
        call MD_outputs%read_time_and_forces_from_files( time_aux, forces )
        call terminate_if_false( time == time_aux, "Last time t is not the same across RT-TDDFT and MD output files")
      end if
    else
      time = real_zero
      e_field%components = real_zero
      p_vec%components = real_zero
    end if
    
    eps_occ = input%groundstate%epsocc

    if( my_rank_writes_to_output ) then
      call open_rttddft_outputs( rt )
      call write_file_info_header()
    end if
    
    call initialize_rttddft( rt, propagator, vec_pot, a_tot_save, molecular_dynamics, &
      psi, overlap, ham_init, ham_time, ham_past, effective_potential_init, &
      apwalm, pmat, pmatmt, rhomt_frozen, rhoir_frozen, occupations, initial_ks_energies, k_dependent_dims, &
      eps_occ, kset, Gkset, Gset, ks_lapwlo_transition_matrix, pws_for_berry_phase, k_ptrs, &
      td_overlap_det, berry_coupling_term, prev_phases, e_field, e_vec_save, j_para_spurious, p_vec_init, energy_gap )
    call check_rttddft_setup( propagator%time_step(), initial_ks_energies, rt%use_berry_phase(), &
      vec_pot, time, rt%t_end, avec, kset%ngridk, energy_gap )
    
    call distribute_loop( mpi_env_k, kset%nkpt, first_kpt, last_kpt )
    dt = propagator%time_step()
    if( molecular_dynamics%on ) then
      call init_MD( rt%do_from_scratch(), time, vec_pot%a_tot, dt, psi%active, &
        occupations, overlap, ham_time, kset, time_step_multiplier, molecular_dynamics, &
        MD_outputs, nuclei_motion, e_field, forces )
      if( rt%restart_previous_calculation() ) then
        call nuclei_motion%allocate_arrays( natmtot )
        call read_state_Ehrenfest_MD( nuclei_motion, rt%restart_file_handler, mpi_env_k )
        call nuclei_motion%update_globals( )
        call update_exciting_globals_for_new_ions_positions( first_kpt, apwalm )
        call initialize_me( Gset )
        if( molecular_dynamics%update_overlap .or. allocated(mathcalH) .or. &
            & allocated(mathcalB) .or. molecular_dynamics%update_pmat ) then
          if( molecular_dynamics%update_pmat ) &
            call obtain_pmat_LAPWLOBasis( first_kpt, rt%pmat%force_pmat_hermitian, apwalm, pmat, pmatmt )
          if( molecular_dynamics%update_overlap ) call overlap%calculate( apwalm, Gkset, &
            pmatmt, vec_pot%a_tot, mathcalH=mathcalH, mathcalB=mathcalB )
          if( propagator%extrapolation_needed() ) ham_past = ham_time  
          call update_hamiltonian_without_pa_term_lapw( first_kpt, lmax_potential, vec_pot%a_tot, &
            ham_time, apwalm, Gkset, update_mathcalH=allocated( mathcalH ) )
          call add_external_coupling_velocity_gauge( vec_pot%a_tot, overlap, ham_time, pmat, k_dependent_dims )
        end if
      end if
    end if

    if( rt%restart_previous_calculation() .and. rt%use_velocity_gauge() ) then
      call j_ind%evaluate_paramagnetic( psi, pmat, occupations, kset%wkpt(first_kpt:last_kpt), mpi_env_k )
      call j_ind%evaluate_diamagnetic( chgval/Omega, vec_pot%a_tot )
    end if

    ! Allocate variables to be stored and printed only after rt_input%n_print steps
    allocate( time_store(rt%n_print), a_ind_store(rt%n_print), a_tot_store(rt%n_print))
    allocate( j_ind_store(rt%n_print), p_vec_store(rt%n_print), e_vec_store(rt%n_print) )
    if( rt%printTimings%general() ) allocate( timing_store(rt%n_print) )
    if( rt%calculate_total_energy ) allocate( etotstore(rt%n_print) )
    if( rt%calculate_n_exc ) allocate( n_exc(rt%n_print), n_gs(rt%n_print) )
    if( molecular_dynamics%on ) then
      allocate( print_forces(rt%n_print) )
      allocate( nuclei_motion_store(rt%n_print) )
      do i_print = 1, rt%n_print
        call nuclei_motion_store(i_print)%allocate_arrays( natmtot )
      end do
      if ( molecular_dynamics%print_all_force_components ) then
        allocate( forces_store(rt%n_print) )
        do i_print = 1, rt%n_print
          call forces_store(i_print)%allocate_arrays( natmtot, .true. )
        end do
      end if
    end if
    e_field%components = - vec_pot%get_dA_dt( time ) / c
    if( my_rank_writes_to_output .and. rt%do_from_scratch() ) &
      call write_fields( [time], [vec_pot%a_ind], [vec_pot%a_tot], [p_vec], [j_ind%total()], [e_field] )

    ! Total energy
    if ( rt%calculate_total_energy .and. rt%do_from_scratch() ) then
      if ( psi%has_frozen() ) call update_density( first_kpt, psi, occupations, 0, &
        .false., rt%l_rad_step, rhomt_frozen, rhoir_frozen, ks_lapwlo_transition_matrix )
      call potcoul()
      call potxc()
      call obtain_energy_rttddft( first_kpt, ham_time, psi, occupations, &
        initial_ks_energies, mpi_env_k, kset%wkpt(first_kpt:last_kpt), etotstore(1) )
      if( my_rank_writes_to_output ) call write_total_energy( .True., [time], [etotstore(1)] )
    end if

    ! Number of excitations
    if ( rt%calculate_n_exc .and. rt%do_from_scratch() ) then
      call psi%obtain_number_excitations( overlap, eps_occ, occupations, kset%wkpt(first_kpt:last_kpt), mpi_env_k, n_exc(1), n_gs(1) )
      if( my_rank_writes_to_output ) call write_nexc( .True., [time], [n_exc(1)], [n_gs(1)] )
    end if

    if ( rt%screenshots%on ) then
      if ( rt%screenshots%density%on ) then
        call update_density( first_kpt, psi, occupations, 0, rt%normalize_WF, rt%l_rad_step, &
          rhomt_frozen, rhoir_frozen, ks_lapwlo_transition_matrix, dens_case = ground_state )
        rhomt_init = rhomt
        rhoir_init = rhoir
      end if
      if( rt%do_from_scratch() ) call screenshot( 0, rt%screenshots, overlap, psi, &
        ham_time, k_dependent_dims, occupations, rhomt, rhoir, mpi_env=mpi_env_k )
    end if ! rt%screenshots%on

    if( rt%printTimings%general() ) then
      call timesec( timef )
      if( my_rank_writes_to_output ) call write_timing( timef - timei ) ! write time for initialization
      timeiter = timef
    end if

    ! whether explicitly field-independent Hamiltonian should be evolved in time
    evolve_H0 = ( molecular_dynamics%on .or. ( .not. rt%eeInteraction%use_ipa() ) )

    i_print = 1
    first_step = int( time / dt, kind = i32 ) + 1
    last_step = int( rt%t_end / dt, kind = i32 )
    ! This is the most important loop (performed for each time step \(\Delta t\)
    do it = first_step, last_step

      call timing%reset()
      ! Variable to store the timing of each iteration
      if( rt%printTimings%general() ) timei = timeiter

      ! Shall the screenshot be taken on the current step
      take_screenshot = .false.
      if ( rt%screenshots%on ) take_screenshot = ( mod( it, rt%screenshots%n_steps ) == 0 ) .or. ( it == last_step )
        
      ! We may need to update charge density on some steps
      density_needed = ( .not. rt%eeInteraction%use_ipa() )
      if ( take_screenshot ) density_needed = density_needed .or. rt%screenshots%density%on

      ! WAVEFUNCTION
      if ( save_wavefunction( rt, i_print, it, last_step, propagator ) ) call psi%save()
      if ( molecular_dynamics%on .and. molecular_dynamics%basis_derivative ) then
        call update_basis_derivative( nuclei_motion%velocities, mathcalB, B_time, B_past )
        ham_time = ham_time - zi*B_time
      end if
      call propagator%evolve( list_of_H_minus_dt=ham_past, list_of_H_0=ham_time, &
        list_of_S=overlap%array, psi=psi%active, dims=k_dependent_dims )
      if ( rt%normalize_WF ) call psi%normalize( overlap )
      if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%wavefunction )

      ! The "real time" t of our evolution is now t_in + dt = t_out, 
      ! where the step is from t_in to t_out
      time = time + dt
      e_field%components = - vec_pot%get_dA_dt( time ) / c

      j_ind_save = j_ind
      p_vec_prev = p_vec
      
      if ( rt%use_velocity_gauge() ) then
        ! Update the paramagnetic component of the induced current density
        call j_ind%evaluate_paramagnetic( psi, pmat, occupations, kset%wkpt(first_kpt:last_kpt), mpi_env_k )
        if ( rt%subtract_J0 ) call j_ind%paramagnetic%add_vector( -j_para_spurious%components )
      else
        call get_td_overlap_det_and_berry_coupling_term( first_kpt, e_field, pws_for_berry_phase, &
          psi, kset, k_ptrs, td_overlap_det, berry_coupling_term )
        call p_vec%get_with_mtp( td_overlap_det, kset%ngridk, kset%ikmap, avec, prev_phases, .true. )
        call p_vec%add_vector( - p_vec_init%components )
        if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%td_berry )
        call p_vec%obtain_j( p_vec_prev, dt, j_ind%paramagnetic )
      end if
      if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%current_density )

      ! DENSITY
      if ( density_needed ) call update_density( first_kpt, psi, occupations, it, rt%normalize_WF, &
        rt%l_rad_step, rhomt_frozen, rhoir_frozen, ks_lapwlo_transition_matrix, rt%printTimings, timing%t_RTTDDFT%dens )
      
      ! KS-POTENTIAL
      if ( .not. rt%eeInteraction%use_ipa() ) call update_potential( rt%printTimings, timing%t_RTTDDFT%pot, rt%eeInteraction%coulomb_only() )
      
      if( rt%printTimings%general() ) call timesec( timei )
      ! Check if we need to save aind, pvec, atot and aext
      if( vec_pot%is_external_field_given() .and. rt%predictor_corrector%on .and. ( .not. vec_pot%is_solver_euler() ) ) then
        a_ind_save = vec_pot%a_ind
        p_vec_save = p_vec
      end if
      a_tot_save = vec_pot%a_tot
      
      call update_a_ind_and_p_vec( time, dt, j_ind_save, j_ind%paramagnetic, vec_pot, p_vec_prev )
      if ( rt%use_velocity_gauge() ) p_vec = p_vec_prev
      call vec_pot%evaluate_a_tot( time )
      if( molecular_dynamics%on ) then
        if( vec_pot%is_total_field_given() ) then
          call e_field%obtain_electric_field( 2*dt, Vector_Potential_Field(vec_pot%applied_vector_potential( time+dt )), a_tot_save )
        else
          call e_field%obtain_electric_field( dt, vec_pot%a_tot, a_tot_save )
        end if
      end if

      if( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%vector_potential )

      if ( rt%use_velocity_gauge() ) call j_ind%evaluate_diamagnetic( chgval/Omega, vec_pot%a_tot )
      if ( rt%predictor_corrector%on .and. ( .not. vec_pot%is_solver_euler() ) ) j_ind_save = j_ind

      ! HAMILTONIAN
      if( propagator%extrapolation_needed() ) ham_past = ham_time
      if ( .not. evolve_H0 ) then
        ham_time = ham_init
      else
        if ( rt%use_lapwlo_basis() ) then
          call update_hamiltonian_without_pa_term_lapw( first_kpt, lmax_potential, vec_pot%a_tot, &
            ham_time, apwalm, Gkset, rt%printTimings, timing%t_RTTDDFT%ham )
        else
          call update_hamiltonian_without_pa_term_ks( first_kpt, lmax_potential, ham_time, &
            apwalm, ks_lapwlo_transition_matrix, effective_potential_init, ham_init, &
            Gkset, rt%printTimings, timing%t_RTTDDFT%ham )
        end if 
      end if
      if ( rt%use_velocity_gauge() ) then
        call add_external_coupling_velocity_gauge( vec_pot%a_tot, overlap, ham_time, pmat, k_dependent_dims )
      else
        call add_external_coupling_berry_phase( berry_coupling_term, ham_time, k_dependent_dims )
      end if

      if ( rt%predictor_corrector%on ) then
        if ( rt%printTimings%general() ) call timesec( timei )
        call loop_predictor_corrector( it, time, rt, first_kpt, psi, occupations, overlap, &
          ham_time, ham_past, k_dependent_dims, apwalm, pmat, a_ind_save, a_tot_save, &
          p_vec_save, j_ind_save, j_para_spurious, propagator, vec_pot, p_vec, j_ind, &
          mpi_env_k, pred_corr_reached_max_steps, lmax_potential, kset, Gkset, &
          pws_for_berry_phase, k_ptrs, td_overlap_det, berry_coupling_term, ks_lapwlo_transition_matrix, &
          effective_potential_init, ham_init, rhomt_frozen, rhoir_frozen )
        if ( pred_corr_reached_max_steps .and. my_rank_writes_to_output ) &
          call warning( 'Problems with convergence (PredCorr), time: ' //  to_char(time) )
        if ( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%pred_corr )
      end if !predictor-corrector

      ! Obtain the total energy, if requested
      if( rt%calculate_total_energy ) then
        if ( rt%printTimings%detailed() ) call timesec( timei )
        call obtain_energy_rttddft( first_kpt, ham_time, psi, occupations, initial_ks_energies, mpi_env_k, &
        kset%wkpt(first_kpt:last_kpt), etotstore(i_print) )
        if ( rt%printTimings%detailed() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%energy )
      end if

      ! Obtain the number of excited electrons, if requested
      if( rt%calculate_n_exc ) then
        if ( rt%printTimings%detailed() ) call timesec( timei )
        call psi%obtain_number_excitations( overlap, eps_occ, occupations, kset%wkpt(first_kpt:last_kpt), mpi_env_k, n_exc(i_print), n_gs(i_print))
        if( rt%printTimings%detailed() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%n_exc )
      end if

      if ( molecular_dynamics%on ) then
        print_forces(i_print) = .False.
        if ( mod( it, time_step_multiplier ) == 0 ) then
          if ( rt%printTimings%general() ) call timesec( timei )
          call forces%save_total_force()
          call force_rttdft( forces, vec_pot%a_tot, e_field, molecular_dynamics, &
            psi%active, occupations, overlap%array, ham_time, &
            kset%wkpt, rt%printTimings, timing%t_Ehrenfest )
          call move_ions( first_kpt, forces%total, forces%total_save, molecular_dynamics%time_step, &
            nuclei_motion, apwalm, rt%printTimings, timing%t_Ehrenfest )
          print_forces(i_print) = .True.
          nuclei_motion_store(i_print) = nuclei_motion
          forces_store(i_print) = forces
          ! Update Hamiltonian with the new basis
          if( molecular_dynamics%update_overlap .or. allocated(mathcalH) .or. &
            & allocated(mathcalB) .or. molecular_dynamics%update_pmat ) then
              if( molecular_dynamics%update_pmat ) then
                if( rt%printTimings%detailed() ) call timesec( time_aux )
                call obtain_pmat_LAPWLOBasis( first_kpt, rt%pmat%force_pmat_hermitian, apwalm, pmat, pmatmt )
                if( rt%printTimings%detailed() ) call timesec_RTTDDFT( time_aux, timing%t_Ehrenfest%pmat )
              end if
            if( propagator%extrapolation_needed() ) ham_past = ham_time  

            call initialize_me( Gset )
            if ( molecular_dynamics%update_overlap ) then
              call initialize_me( Gset )
              call overlap%calculate( apwalm, Gkset, pmatmt, vec_pot%a_tot, timing%t_Ehrenfest%overlap, &
                mathcalH=mathcalH, mathcalB=mathcalB )
            end if

            call update_hamiltonian_without_pa_term_lapw( first_kpt, lmax_potential, vec_pot%a_tot, ham_time, apwalm, &
              Gkset,rt%printTimings, timing%t_RTTDDFT%ham, timing%t_Ehrenfest, update_mathcalH=allocated( mathcalH ) )
            call add_external_coupling_velocity_gauge( vec_pot%a_tot, overlap, ham_time, pmat, k_dependent_dims )
          end if
          if( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_Ehrenfest%t_MD_step )
        end if ! if ( mod( it, time_step_multiplier ) == 0 )
      end if ! if ( molecular_dynamics%on ) then

      ! Check if a screenshot has been requested
      if ( take_screenshot ) then
        if( rt%printTimings%general() ) call timesec( timei )
        call screenshot( it, rt%screenshots, overlap, psi, ham_time, k_dependent_dims, &
          occupations, rhomt, rhoir, rhomt_init, rhoir_init, mpi_env_k )
        if( rt%printTimings%general() ) call timesec_RTTDDFT( timei, timing%t_RTTDDFT%screenshot )
      end if

      ! Store relevant information from this iteration
      time_store(i_print) = time
      a_ind_store(i_print) = vec_pot%a_ind
      a_tot_store(i_print) = vec_pot%a_tot
      j_ind_store(i_print) = j_ind%total()
      p_vec_store(i_print) = p_vec
      e_vec_store(i_print) = e_field
      if( rt%printTimings%general() ) timing_store(i_print) = timing

      ! Print relevant information, every 'rt%n_print' steps
      if ( i_print == rt%n_print ) then
        call timesec( timei )
        if( my_rank_writes_to_output ) then
          call write_fields( time_store, a_ind_store, a_tot_store, p_vec_store, j_ind_store, e_vec_store )
          if ( rt%calculate_total_energy ) call write_total_energy( .False., time_store, etotstore )
          if ( rt%calculate_n_exc ) call write_nexc( .False., time_store, n_exc, n_gs )

          if( molecular_dynamics%on ) then
            do i_print = 1, rt%n_print
              if( print_forces(i_print) ) call MD_outputs%write_to_files( &
                time_store(i_print), nuclei_motion_store(i_print), forces_store(i_print) )
            end do
          end if
        end if
        if( rt%write_restart() ) then
          call write_wavefunction( t, first_kpt, kset%vkl(:, first_kpt:last_kpt), &
            psi%active, mpi_env_k, rt%restart_file_handler, kset%nkpt )
          if( propagator%extrapolation_needed() ) &
            call write_wavefunction( t_minus_dt, first_kpt, kset%vkl(:, first_kpt:last_kpt), &
              psi%active_save, mpi_env_k, rt%restart_file_handler, kset%nkpt )
          if( molecular_dynamics%on ) call write_state_Ehrenfest_MD( nuclei_motion, rt%restart_file_handler, mpi_env_k )
          if ( rt%use_berry_phase() ) then
            if ( my_rank_writes_to_output ) call write_phases( prev_phases, rt%restart_file_handler, mpi_env_k )
          end if
        end if
        if( rt%printTimings%general() ) then
          call timesec_RTTDDFT( timei, timing_store(rt%n_print)%t_RTTDDFT%t_print )
          call timesec_RTTDDFT( timeiter, timing_store(rt%n_print)%t_iteration )
          if( my_rank_writes_to_output ) call write_timing( it, timing_store, molecular_dynamics%on )
        end if
        i_print = 0
      else ! if ( iprint == rt%n_print ) 
        if( rt%printTimings%general() ) call timesec_RTTDDFT( timeiter, timing_store(i_print)%t_iteration )
      end if ! if ( iprint == rt%n_print ) 
      i_print = i_print + 1
      ! Make all the processes wait here: the master alone has been writing the files above
      call barrier( mpi_env_k )
    end do ! it = first_step, last_step

    ! write outputs, if i_print is not 1 (which means that last_step is not a multiple of rt%n_print)
    if( i_print /= 1 ) then
      if( my_rank_writes_to_output ) then
        associate( n => i_print-1 )
          call write_fields( time_store(1:n), a_ind_store(1:n), a_tot_store(1:n), p_vec_store(1:n), j_ind_store(1:n), e_vec_store(1:n) )
          if ( rt%calculate_total_energy ) call write_total_energy( .False., time_store(1:n), etotstore(1:n) )
          if ( rt%calculate_n_exc ) call write_nexc( .False., time_store(1:n), n_exc(1:n), n_gs(1:n) )
          if( molecular_dynamics%on ) then
            do it = 1, n
              if( print_forces(it) ) call MD_outputs%write_to_files( time_store(it), nuclei_motion_store(it), forces_store(it) )
            end do
          end if
        end associate
      end if
      call write_wavefunction( t, first_kpt, kset%vkl(:, first_kpt:last_kpt), psi%active, &
        mpi_env_k, rt%restart_file_handler, kset%nkpt )
      if( propagator%extrapolation_needed() ) call write_wavefunction( t_minus_dt, first_kpt, &
        kset%vkl(:, first_kpt:last_kpt), psi%active_save, mpi_env_k, rt%restart_file_handler, kset%nkpt )
      if ( rt%use_berry_phase() ) then
        if ( my_rank_writes_to_output ) call write_phases( prev_phases, rt%restart_file_handler, mpi_env_k )
      end if
    end if
    if ( molecular_dynamics%on ) then 
      call write_state_Ehrenfest_MD( nuclei_motion, rt%restart_file_handler, mpi_env_k )
      call deallocate_global_arrays()
    end if
    call me_finit()

    if ( my_rank_writes_to_output ) then
      call write_file_info( 'Real-time TDDFT calculation finished' )
      call close_rttddft_outputs( rt )
      if ( molecular_dynamics%on ) call MD_outputs%close_files()
    end if
    
  end subroutine coordinate_rttddft_calculation

  !> (private) Check if we should save the current wavefunction coefficients in a new array
  pure logical function save_wavefunction( rt_inp, i_print, it, last_step, prop ) result( check )
    !> Type that encapsulates the input keywords
    type(rttddft_input_keys), intent(in) :: rt_inp
    !> Counter for printing out outputs
    integer(i32), intent(in) :: i_print
    !> Counter for RT-TDDFT iteration steps
    integer(i32), intent(in) :: it
    !> Index of the last step in RT-TDDFT iterations
    integer(i32), intent(in) :: last_step
    !> Propagator used
    class(propagator_type), intent(in) :: prop

    logical :: is_print_step, is_last_step, print_wavefunction_now

    is_print_step = ( i_print == rt_inp%n_print )
    is_last_step = ( it == last_step )
    print_wavefunction_now = ( rt_inp%write_restart() .and. is_print_step ) .or. is_last_step
    check = prop%extrapolation_needed() .and. print_wavefunction_now
  end function

  !> Open output files
  subroutine open_rttddft_outputs( rt_input )
    !> Type that encapsulates the input keywords
    type(rttddft_input_keys), intent(in) :: rt_input

    call open_files_vector_fields( new=rt_input%do_from_scratch() )
    call open_file_info()
    if( rt_input%calculate_total_energy ) call open_file_etot( new=rt_input%do_from_scratch() )
    if( rt_input%calculate_n_exc ) call open_file_nexc( new=rt_input%do_from_scratch() )
    if( rt_input%printTimings%general() ) call open_file_timing( new=rt_input%do_from_scratch() )
  end subroutine

  !> Close output files
  subroutine close_rttddft_outputs( rt_input )
    !> Type that encapsulates the input keywords
    type(rttddft_input_keys), intent(in) :: rt_input

    call close_files_vector_fields()
    call close_file_info()
    if( rt_input%calculate_total_energy ) call close_file_etot()
    if( rt_input%calculate_n_exc ) call close_file_nexc()
    if( rt_input%printTimings%general() ) call close_file_timing()
  end subroutine

  !> (private) Read time and fields stored in the corresponding files
  subroutine read_time_and_fields( t, p_vec_t, a_t, a_ind_t_minus_dt, a_tot_t_minus_dt, &
      e_vec_t, e_vec_t_minus_dt )
    !> Time \(t\)
    real(dp), intent(out) :: t
    !> Polarization vector at time \(t\)
    type(Polarization), intent(out) :: p_vec_t
    !> Vector potential at time \(t\)
    type(Vector_Potential), intent(inout) :: a_t
    !> \(\mathbf{A}_{ind}) at time \(t-\Delta t\)
    type(Vector_Potential_Field), intent(out) :: a_ind_t_minus_dt
    !> \(\mathbf{A}_{tot}) at time \(t-\Delta t\)
    type(Vector_Potential_Field), intent(out) :: a_tot_t_minus_dt
    !> Electric field at time \(t\)
    type(Electric_Field), intent(out) :: e_vec_t
    !> Electric field at time \(t-\Delta t\)
    type(Electric_Field), intent(out) :: e_vec_t_minus_dt

    call read_vector_field( t, p_vec_t )
    call read_vector_field( t, e_vec_t, e_vec_t_minus_dt )
    call read_vector_field( t, a_t%a_ind, a_ind_t_minus_dt, a_t%a_tot, a_tot_t_minus_dt )
  end subroutine

  !> Wrapper to call [[write_vector_field]]
  subroutine write_fields( time_array, a_ind_array, a_tot_array, p_vec_array, j_ind_array, e_vec_array )
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
    !> Array with the external electric field
    type(Electric_Field), intent(in) :: e_vec_array(:)

    call write_vector_field( time_array, a_ind_array, a_tot_array )
    call write_vector_field( time_array, p_vec_array )
    call write_vector_field( time_array, j_ind_array )
    call write_vector_field( time_array, e_vec_array )
  end subroutine

  
  !> Loop used in the predictor-corrector method
  subroutine loop_predictor_corrector( it, time, rt, first_kpt, psi, occupations, &
    overlap, ham_time, ham_past, k_dependent_dims, apwalm, pmat, &
    a_ind_t_minus_dt, a_tot_t_minus_dt, p_vec_t_minus_dt, j_t_minus_dt, j_para_spurious,&
    propagator, a_t, p_vec, j_t, mpi_env, max_steps_reached, lmax_potential, kset, Gkset, &
    pws_for_berry_phase, k_ptrs, td_overlap_det, berry_coupling_term, ks_lapwlo_transition_matrix, &
    effective_potential_init, ham_init, rhomt_frozen, rhoir_frozen )
    !> current iteration number in the RT-TDDFT loop
    integer(i32), intent(in) :: it
    !> time \( t \)
    real(dp), intent(in) :: time
    !> Type that encapsulates the parameters defined in the input file
    type(rttddft_input_keys), intent(in) :: rt
    !> index of the first `k-point` to be considered in the sum
    integer(i32), intent(in) :: first_kpt
    !> Basis-expansion coefficients of the KS-WFs
    class(wavefunction_set), intent(inout) :: psi
    !> Initial occupations array
    real(dp), intent(in) :: occupations(:, :)
    !> Overlap matrix (of basis functions)
    class(overlap_set), intent(in) :: overlap
    !> Hamiltonian matrix at current time \(t\)
    complex(dp), contiguous, intent(inout) :: ham_time(:, :, :)
    !> Hamiltonian matrix at previous time \(t - \Delta t \)
    complex(dp), contiguous, intent(in) :: ham_past(:, :, :)
    !> k-dependent Hamiltonian dimensions array
    integer(i32), contiguous, intent(in) :: k_dependent_dims(:)
    !> Matching coefficients of the (L)APWs
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
    complex(dp), contiguous, intent(inout) :: pmat(:, :, :, :)
    !> `a_ind` at time \( t-\Delta t\) 
    class(Vector_Potential_Field), intent(in) :: a_ind_t_minus_dt
    !> `a_tot` at time \( t-\Delta t\) 
    class(Vector_Potential_Field), intent(in) :: a_tot_t_minus_dt
    !> Polarization at time \( t-\Delta t\) 
    type(Polarization), intent(in) :: p_vec_t_minus_dt
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
    !> Current density at time \( t \) 
    type(Current_Density), intent(inout) :: j_t
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> When `.True.`, it informs that the maximum steps have been reached
    logical, intent(out) :: max_steps_reached
    !> Maximum value of l used for the potential expansion in MT
    integer(i32), intent(in) :: lmax_potential
    !> Set of k vectors used in the module
    type(k_set), intent(in) :: kset
    !> Set of G+k vectors used for the matrix elements evaluation
    type(Gk_set), intent(in) :: Gkset
    !> Planewave matrix elements between neighbouring \( \mathbf{k} \) points
    complex(dp), contiguous, intent(in) :: pws_for_berry_phase(:, :, :, :, :)
    !> Array containing indices of the neighbouring \( \mathbf{k} \) points
    integer(i32), contiguous, intent(in) :: k_ptrs(:, :, :)
    !> Determinants of the time-dependent overlaps of the periodic parts of the KS-Bloch states
    complex(dp), contiguous, intent(out) :: td_overlap_det(:, :)
    !> Field coupling with the external field constructed with dynamical Berry phase approach
    complex(dp), contiguous, intent(out) :: berry_coupling_term(:, :, :)
    ! KS-LAPW+lo transition matrix (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: ks_lapwlo_transition_matrix(:, :, :)
    !> Effective potential matrix at time \(t = 0 \) 
    complex(dp), contiguous, optional, intent(in) :: effective_potential_init(:, :, :)
    !> Hamiltonian matrix at time \(t = 0 \)
    complex(dp), contiguous, optional, intent(in) :: ham_init(:, :, :)
    !> Frozen part of the muffin-tin density (lmmaxvr, nrmtmax, natmtot)
    real(dp), contiguous, optional, intent(in) :: rhomt_frozen(:, :, :)
    !> Frozen part of the IR density (ngrtot)
    real(dp), contiguous, optional, intent(in) :: rhoir_frozen(:)

    integer(i32) :: i, last_kpt
    real(dp) :: err, dt
    complex(dp), allocatable :: ham_predcorr(:, :, :)
    type(Electric_Field) :: e_vec

    dt = propagator%time_step()
    last_kpt = first_kpt + size( ham_time, 3 ) - 1
    allocate( ham_predcorr, mold = ham_time )

    do i = 1, rt%predictor_corrector%max_steps
      ! WAVEFUNCTION
      call psi%restore()
      call propagator%evolve( list_of_H_dt=ham_time, list_of_H_0=ham_past, list_of_S=overlap%array, psi=psi%active, dims=k_dependent_dims )
      if ( rt%normalize_WF ) call psi%normalize( overlap )

      ! Update the paramagnetic component of the induced current density
      j_t = j_t_minus_dt
      if ( rt%use_velocity_gauge() ) then
        call j_t%evaluate_paramagnetic( psi, pmat, occupations, kset%wkpt(first_kpt:last_kpt), mpi_env )
        if ( rt%subtract_J0 ) call j_t%paramagnetic%add_vector( -j_para_spurious%components )
      else
        e_vec%components = - a_t%get_dA_dt( time ) / c
        call get_td_overlap_det_and_berry_coupling_term( first_kpt, e_vec, pws_for_berry_phase, &
          psi, kset, k_ptrs, td_overlap_det, berry_coupling_term )
      end if

      ! DENSITY
      call update_density( first_kpt, psi, occupations, it, rt%normalize_WF, rt%l_rad_step, &
        rhomt_frozen, rhoir_frozen, ks_lapwlo_transition_matrix )
      ! KS-POTENTIAL
      call update_potential( coulomb_only = rt%eeInteraction%coulomb_only() )

      if ( rt%use_velocity_gauge() ) then
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
      end if

      ! HAMILTONIAN
      ham_predcorr = ham_time
      if ( rt%use_lapwlo_basis() ) then
        call update_hamiltonian_without_pa_term_lapw( first_kpt, lmax_potential, a_t%a_tot, ham_time, apwalm, Gkset )
      else
        call update_hamiltonian_without_pa_term_ks( first_kpt, lmax_potential, ham_time, &
          apwalm, ks_lapwlo_transition_matrix, effective_potential_init, ham_init, Gkset )
      end if
      if ( rt%use_velocity_gauge() ) then
        call add_external_coupling_velocity_gauge( a_t%a_tot, overlap, ham_time, pmat, k_dependent_dims )
      else
        call add_external_coupling_berry_phase( berry_coupling_term, ham_time, k_dependent_dims )
      end if

      ! Check the difference between the two hamiltonians
      err = maxval( abs( ham_predcorr - ham_time ) )
      if ( err <= rt%predictor_corrector%tol ) exit

    end do
    max_steps_reached = ( i > rt%predictor_corrector%max_steps )
  end subroutine 

  !> Subroutine to initialize all MD related variables
  subroutine init_MD( from_scratch, t_0, a_tot, timeStepRTTDDFT, evecfv_time, occupations, overlap, ham_time, &
    kset, time_step_multiplier, molecular_dynamics, MD_outputs, nuclei_motion, e_vec, forces )
    !> If `.true.`, this calculation is done from scratch (i.e. it is not restarting a previous calculation)
    logical, intent(in) :: from_scratch
    !> Initial time \( t_0 \)
    real(dp), intent(in) :: t_0
    !> Vector potential (total)
    class(Vector_Potential_Field), intent(in) :: a_tot
    !> Time step used in the real-time TDDFT calculation
    real(dp), intent(in) :: timeStepRTTDDFT
    !> Basis-expansion coefficients of the KS-WFs at time \(t\)
    complex(dp), intent(in) :: evecfv_time(:, :, :)
    !> State occupations array
    real(dp), intent(in) :: occupations(:, :)
    !> Overlap matrix (of basis functions)
    class(overlap_set), intent(in) :: overlap
    !> Hamiltonian matrix at current time \(t\)
    complex(dp), intent(in) :: ham_time(:, :, :)
    !> k set used in the RT TDDFT module
    type(k_set), intent(in) :: kset
    !> Integer ratio between the time step used in MD and `timeStepRTTDDFT`
    integer(i32), intent(out) :: time_step_multiplier
    !> variable with interfaces to elements defined in the input file
    type(MD_input_keys), intent(inout) :: molecular_dynamics
    !> variable with interfaces to MD outputs
    type(MD_out), intent(inout) :: MD_outputs
    !> This argument packs nuclei positions and velocities
    class(trajectory), intent(inout) :: nuclei_motion
    !> Electric field
    type(Electric_Field), intent(in) :: e_vec
    !> forces acting on all atoms
    type(force), intent(inout) :: forces


    time_step_multiplier = int( molecular_dynamics%time_step/timeStepRTTDDFT )
    molecular_dynamics%time_step = time_step_multiplier*timeStepRTTDDFT
    
    call MD_allocate_global_arrays( nspecies )
    call MD_evaluate_charge_val( )
    
    call forces%allocate_arrays( natmtot, from_scratch )
    if( from_scratch ) then
      call force_rttdft( forces, a_tot, e_vec, molecular_dynamics, evecfv_time, &
        occupations, overlap%array, ham_time, kset%wkpt )
      call nuclei_motion%allocate_arrays( natmtot )
      call nuclei_motion%initialize( input%structure )
      if( molecular_dynamics%basis_derivative ) call update_basis_derivative( nuclei_motion%velocities, mathcalB, B_time, B_past )
    end if
    
    if ( mpiglobal%is_root ) then
      call MD_outputs%open_files( from_scratch, natmtot, molecular_dynamics%print_all_force_components  )
      if( from_scratch ) call MD_outputs%write_to_files( t_0, nuclei_motion, forces )
    end if

  end subroutine

  !> Deallocate global arrays
  subroutine deallocate_global_arrays( )
    call MD_deallocate_global_arrays
    if ( allocated( mathcalH ) ) deallocate( mathcalH )
    if ( allocated( mathcalB ) ) deallocate( mathcalB )
    if ( allocated( B_past ) ) deallocate( B_past )
    if ( allocated( B_time ) ) deallocate( B_time )
  end subroutine
end module rttddft_main
