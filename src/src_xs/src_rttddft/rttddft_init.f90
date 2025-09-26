! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Created Jan 2021 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module implementing general initializations for RT-TDDFT
module rttddft_init
  use asserts, only: assert
  use constants, only: zzero, real_zero, zone, zi
  use exciting_mpi, only: xmpi_bcast, xmpi_allreduce, xmpi_allgatherv
  use general_find_vbm_cbm, only: find_vbm_cbm
  use m_gndstateq, only: gndstateq
  use math_utils, only: plane_wave_in_spherical_harmonics
  use matrix_elements, only: me_init, me_mt_prepare, me_ir_alloc, me_ir_prepare, &
    me_mt_mat, me_ir_mat, me_mt_alloc
  use MD, only: MD_input_keys
  use mod_APW_LO, only: apwordmax, apword, nlorb, lorbl, lofr, apwfr
  use mod_atoms, only: natmtot, spr, nspecies, idxas, atposc, natoms
  use mod_bands, only: evalfv, nomax, numin, ikcbm, ikvbm, ikvcm
  use mod_core_states, only: ncg
  use mod_corestate, only: rhocr, evalcr
  use mod_eigensystem, only: nmatmax, nmat
  use mod_eigenvalue_occupancy, only: nstfv, efermi
  use mod_gvector, only: intgv, ngvec, sfacg, vgc
  use mod_gkvector, only: ngk, ngkmax, gkc, tpgkc, sfacgk, gkmax, vgkc
  use mod_kpointset, only: G_set, generate_G_vectors, k_set, &
    generate_k_vectors, Gk_set, generate_Gk_vectors
  use mod_lattice, only: bvec, avec
  use mod_misc, only: filext
  use mod_muffin_tin, only: lmmaxapw, nrmt
  use mod_potential_and_density, only: rhomt, rhoir
  use mod_spin, only: nspnfv
  use modinput, only: input, getstructHybrid, emptynode, input_type
  use modmpi, only: rank, mpi_env_k, barrier, distribute_loop, terminate_if_false, &
    procofindex, firstk, lastk
  use modxs, only: isreadstate0
  use muffin_tin_basis, only: mt_basis_type
  use precision, only: dp, i32, str_128, sp
  use propagators, only: propagator_type => propagator, create_propagator
  use rttddft_berry, only: get_td_overlap_det_and_berry_coupling_term
  use rttddft_CurrentDensity, only: Current_Density, Current_Density_Field
  use rttddft_Density, only: update_density, save_and_frozen, frozen, ground_state
  use rttddft_electric_field, only: Electric_Field
  use rttddft_file_names, only: filename_avec, filename_evec, filename_jind, filename_pvec, RTTDDFT_GND_sufix
  use rttddft_GlobalMDVariables, only: B_past, B_time, mathcalB
  use rttddft_Hamiltonian, only: hamiltonian_set
  use rttddft_hybrids, only: hybrids_used, Set_Dimension_mixed_product_basis, set_barecoul_basis
  use rttddft_input, only: rttddft_input_keys
  use rttddft_io, only: file_pmat_exists, read_pmat, write_pmat, file_pmat_mt_exists, &
    read_pmat_mt, write_pmat_mt, write_file_info, write_file_info_fill_line_with_char, &
    get_filename_pmat, get_filename_pmat_mt, read_wavefunction, groundstate, t, t_minus_dt, &
    read_phases
  use rttddft_Overlap, only: overlap_set
  use rttddft_pmat, only: obtain_pmat_LAPWloBasis, obtain_pmat_KSBasis
  use rttddft_Polarization, only: Polarization
  use rttddft_potential, only: update_potential
  use rttddft_VectorPotential, only: Vector_Potential, Vector_Potential_Field
  use rttddft_Wavefunction, only: wavefunction_set, initialize_wavefunction_set
  use to_char_conversion, only: to_char
  use mod_gen_lo, only: genlofr
  use m_zfftifc, only: zfftifc

  implicit none
  
  private
  
  public :: initialize_rttddft, initialize_me

contains
!> This subroutine initializes many variables used in a RT-TDDFT calculation.
subroutine initialize_rttddft( rt_inp, propagator, vec_pot, a_tot_t_minus_dt, molecular_dynamics, &
    psi, overlap, H, apwalm, pmat, pmatmt, rhomt_frozen, rhoir_frozen, occupations, &
    occs_tol, kset_rttddft, Gkset, Gset, psi_gnd_lapwlo, pws_for_berry_phase, k_ptrs, &
    td_overlap_det, berry_coupling_term, prev_phases, e_vec, e_vec_save, j_para_spurious, p_vec_init, &
    energy_gap )
  !> Argument that encapsulates the input options of rttddft
  type(rttddft_input_keys), intent(in) :: rt_inp
  !> Argument that encapsulates the propagator
  class(propagator_type), allocatable, intent(out) :: propagator
  !> Argument that encapsulates the vector potential
  type(Vector_Potential), intent(in) :: vec_pot
  !> \(\mathbf{A}_{tot}\) at time \( t-\Delta t\) 
  class(Vector_Potential_Field), intent(in) :: a_tot_t_minus_dt
  !> variable that is an interface to the input keys defined in `input.xml` inside the `MD` block
  type(MD_input_keys), intent(in) :: molecular_dynamics
  !> Basis-expansion coefficients of the KS-WFs
  class(wavefunction_set), allocatable, intent(out) :: psi
  !> Type that encapsulates the overlap matrix (of basis functions, to be allocated in `array_allocation` block)
  type(overlap_set), intent(inout) :: overlap
  !> Type that encapsulates the hamiltonian matrices (to be allocated in `array_allocation` block)
  type(hamiltonian_set), intent(out) :: H
  !> Matching coefficients of the (L)APWs (to be allocated in `array_allocation` block)
  complex(dp), allocatable, intent(out) :: apwalm(:, :, :, :, :)
  !> Momentum matrix elements (to be allocated in `array_allocation` block)
  complex(dp), allocatable, intent(out) :: pmat(:, :, :, :)
  !> Muffin-tin part of the momentum matrix (to be allocated in `array_allocation` block)
  complex(dp), allocatable, intent(out) :: pmatmt(:, :, :, :, :)
  !> Frozen part of the muffin-tin density (to be allocated in `array_allocation` block)
  real(dp), allocatable, intent(out) :: rhomt_frozen(:, :, :)
  !> Frozen part of the IR density (to be allocated in `array_allocation` block)
  real(dp), allocatable, intent(out) :: rhoir_frozen(:)
  !> State occupations array (to be allocated in `array_allocation` block)
  real(dp), allocatable, intent(out) :: occupations(:, :)
  !> Minimal value of occupation for the state to be 'occupied'
  real(dp), intent(in) :: occs_tol
  !> Set of \( \mathbf{k} \) points
  type(k_set), intent(out) :: kset_rttddft
  !> Set of \( \mathbf{G} + \mathbf{k} \) vectors used for the matrix elements evaluation
  type(Gk_set), intent(out) :: Gkset
  !> Set of \( \mathbf{G} \) vectors used for the matrix elements evaluation
  type(G_set), intent(out) :: Gset
  !> KS-LAPW+lo transition matrix (ground state set in the LAPW+lo basis) (to be allocated in `array_allocation` block)
  complex(dp), allocatable, intent(out) :: psi_gnd_lapwlo(:, :, :)
  !> Planewave matrix elements between neighbouring \( \mathbf{k} \) points (to be allocated in `array_allocation` block)
  complex(dp), allocatable, intent(out) :: pws_for_berry_phase(:, :, :, :, :)
  !> Array containing indices of the neighbouring \( \mathbf{k} \) points (to be allocated in `array_allocation` block)
  integer(i32), allocatable, intent(out) :: k_ptrs(:, :, :)
  !> Determinants of the time-dependent overlaps of the periodic parts of the KS-Bloch states (to be allocated in `array_allocation` block)
  complex(dp), allocatable, intent(out) :: td_overlap_det(:, :)
  !> Field coupling with the external field constructed with dynamical Berry phase approach
  complex(dp), allocatable, intent(out) :: berry_coupling_term(:, :, :)
  !> When using dynamical Berry phase approach, phases used for the MTP polarization evaluation
  real(dp), allocatable, intent(out) :: prev_phases(:, :)
  !> \(\mathbf{E}\) at time \( t = t_{\rm start} \) 
  type(Electric_Field), intent(in) :: e_vec
  !> \(\mathbf{E}\) at time \( t = t_{\rm start} - \Delta t \) 
  type(Electric_Field), intent(in) :: e_vec_save
  !> Spurious paramagnetic current density obtained at \( t = 0 \)
  type(Current_Density_Field), intent(out) :: j_para_spurious
  !> GS polarization obtained at \( t = 0 \)
  type(Polarization), intent(out) :: p_vec_init
  !> Energy gap
  real(dp), intent(out) :: energy_gap

  integer(i32) :: ik, first_kpt, last_kpt, l_max_pot, ham_dimension, kgrid_neighbours
  logical :: allocate_H0, evolve_H0, my_rank_writes_to_output, success
  type(Vector_Potential_Field) :: a_aux
  type(Current_Density) :: j_aux
  type(Electric_Field) :: e_aux
  complex(dp), allocatable :: psi_gnd_lapwlo_copy(:, :, :)
  integer(i32), allocatable :: ik_to_array_position(:), k_shifts(:, :, :), shift_positions(:)
  real(dp), allocatable :: dk_vec(:, :, :)
  logical, allocatable :: k_needed(:), proc_needed(:)

  call adjust_input_and_init_exciting_globals( input )
  l_max_pot = input%groundstate%lmaxvr

  my_rank_writes_to_output = (rank == 0)
  if ( my_rank_writes_to_output ) then
    call write_file_info_fill_line_with_char('=')
    call write_file_info('Non-self-consistent GS for TDDFT calculations - started'//new_line( 'a' ))
  end if
  
  ! Read from STATE.OUT exclusively
  isreadstate0 = .true.
  ! One-shot GS calculation. Since an XS calculation with Hybrid functionals uses the GS parameters, 
  ! a one shot GS calculation serves no purpose
  if ( .not. hybrids_used() ) call gndstateq( input%xs%vkloff, RTTDDFT_GND_sufix//filext )

  ! Generate k, G, G+k vectors
  call generate_k_vectors( kset_rttddft, bvec, input%groundstate%ngridk, input%xs%vkloff, .false., .false. )
  call generate_G_vectors( Gset, bvec, intgv, input%groundstate%gmaxvr )
  call generate_Gk_vectors( Gkset, kset_rttddft, Gset, gkmax )
  ! Distribute k-points over MPI ranks
  call distribute_loop( mpi_env_k, kset_rttddft%nkpt, first_kpt, last_kpt )
  
  evolve_H0 = ( molecular_dynamics%on .or. ( .not. rt_inp%eeInteraction%use_ipa() ) )
  allocate_H0 = (.not. evolve_H0) .and. rt_inp%use_lapwlo_basis()

  ! Neighbour is a k point from another MPI rank, which is reachable in 1 or 2 jumps
  ! by one of the k points controlled by the current MPI rank
  kgrid_neighbours = 0
  if ( rt_inp%use_berry_phase() ) call get_kgrid_neighbours_info( first_kpt, last_kpt, &
    kset_rttddft, Gset, k_ptrs, dk_vec, k_needed, proc_needed, kgrid_neighbours, &
    shift_positions, k_shifts )

  ham_dimension = nmatmax
  if ( rt_inp%use_ks_basis() ) ham_dimension = nstfv

  call create_propagator( propagator, rt_inp%propagator_input, .not. molecular_dynamics%on )
  array_allocation: block
    allocate( psi_gnd_lapwlo(nmatmax, nstfv, first_kpt : last_kpt + kgrid_neighbours), source = zzero )
    allocate( occupations(nstfv, first_kpt : last_kpt), source = real_zero )
    allocate( apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt + kgrid_neighbours) )
    call H%allocate( ham_dimension, first_kpt, nmat(1, first_kpt:last_kpt), nstfv, &
       propagator%extrapolation_needed(), allocate_H0, rt_inp%use_lapwlo_basis(), molecular_dynamics%valence_corrections )
    call overlap%allocate( rt_inp%use_lapwlo_basis(), ham_dimension, first_kpt, last_kpt )
    if ( rt_inp%use_velocity_gauge() ) allocate( pmat(ham_dimension, ham_dimension, 3, first_kpt : last_kpt) )
    if ( ( molecular_dynamics%on ) .and. ( molecular_dynamics%valence_corrections .or. molecular_dynamics%basis_derivative ) ) &
      allocate( pmatmt(nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt) )
    if ( rt_inp%n_frozen > 0 ) then
      allocate( rhomt_frozen, mold = rhomt )
      allocate( rhoir_frozen, mold = rhoir )
    end if
    if ( rt_inp%use_berry_phase() ) then
      allocate( pws_for_berry_phase(nstfv, nstfv, first_kpt : last_kpt, 3, 4), source = zzero )
      allocate( td_overlap_det(3, kset_rttddft%nkpt), source = zzero )
      allocate( berry_coupling_term, mold = H%H_t%array )
      allocate( prev_phases(maxval( kset_rttddft%ngridk )**2, 3), source = real_zero )
    end if
  end block array_allocation

  call read_WF_potential_rttddft( first_kpt, kset_rttddft, psi_gnd_lapwlo(:, :, first_kpt : last_kpt), &
    occupations, H%initial_eigenvalues )
  ! Attention: `me_init` must be called after `read_WF_potential_rttddft`, since it initialize radial functions
  call initialize_me( Gset )
  
  do ik = first_kpt, last_kpt
    ! Matching coefficients (apwalm)
    call match( ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), sfacgk(:, :, 1, ik), apwalm(:, :, :, :, ik) )
  end do
  ! get apwalm and initiall states from processes with the neighbouring \( \mathbf{k} \) points
  if ( rt_inp%use_berry_phase() ) then
    call get_data_from_neighbours( first_kpt, last_kpt, kset_rttddft%nkpt, &
      k_needed, apwalm, psi_gnd_lapwlo, ik_to_array_position )

    call calc_planewave_matrix_elements( first_kpt, k_ptrs, dk_vec, apwalm, psi_gnd_lapwlo, &
      pws_for_berry_phase, ik_to_array_position, Gkset, Gset, shift_positions, k_shifts )
    
    if ( kgrid_neighbours > 0 ) then
      ! extended versions of apwalm and psi_gnd_lapwlo are not needed anymore
      deallocate( apwalm )
      allocate( psi_gnd_lapwlo_copy, source = psi_gnd_lapwlo(:, :, first_kpt : last_kpt) )
      deallocate( psi_gnd_lapwlo )
      allocate( psi_gnd_lapwlo, source = psi_gnd_lapwlo_copy )
      deallocate( psi_gnd_lapwlo_copy )
      ! apwalm is a heavy array, it is cheaper to rematch it
      allocate( apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt) )
      do ik = first_kpt, last_kpt
        call match(ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), sfacgk(:, :, 1, ik), apwalm(:, :, :, :, ik))
      end do
    end if
  end if ! dynamical Berry phase approach

  call initialize_wavefunction_set( psi, rt_inp%use_lapwlo_basis(), &
    propagator%extrapolation_needed() .or. rt_inp%restart_previous_calculation() , &
    rt_inp%n_frozen, psi_gnd_lapwlo, occupations, occs_tol )
  if ( rt_inp%use_lapwlo_basis() ) deallocate( psi_gnd_lapwlo )

  ! In general, non-physical parameters (such as the k-grid) can differ between the GS
  ! and RT modules, which can result in e.g. different XC potential calculated from 
  ! the same electron density. For consistency, we generate the initial density and potential
  ! at step 0 the same way as during the time propagation. 
  call update_density( first_kpt, psi, occupations, -1, rt_inp%normalize_WF, &
    rt_inp%l_rad_step, rhomt_frozen, rhoir_frozen, psi_gnd_lapwlo, dens_case = ground_state )
  call update_potential()

  ! A special case of an input parameter for the EH and EHM propagators:
  ! first, n_eigvecs_houston can be < 0 and should be redefined as soon as nstfv is known
  ! second, n_eigvecs_houston should be at least n_occupied, and not larger than basis size
  ! The routine does nothing if different propagator is used
  call propagator%update_and_check( psi%n_occupied(), minval( nmat(1, first_kpt : last_kpt) ), nstfv, success )
  call terminate_if_false( success, &
    'Error: Provided value of nEigenvectorsEH is either smaller than the number of occupied states or larger than basis size.' )

  call get_energy_gap( first_kpt, H%initial_eigenvalues, occupations, kset_rttddft%nkpt, energy_gap )

  if( molecular_dynamics%on ) call allocate_MD_globals( first_kpt, last_kpt, &
    allocate_B=molecular_dynamics%basis_derivative, &
    allocate_mathcalB=molecular_dynamics%valence_corrections .or. molecular_dynamics%basis_derivative )  
  
  if ( my_rank_writes_to_output ) call estimate_memory_and_write_to_info( molecular_dynamics%on, &
    rt_inp%predictor_corrector%on, psi, overlap, H, apwalm, kgrid_neighbours, &
    kset_rttddft%nkpt, pmat, pmatmt, psi_gnd_lapwlo, pws_for_berry_phase, berry_coupling_term )

  if ( hybrids_used() ) then
    if ( input%xs%realTimeTDDFT%calcNonlocalCurrentDensity ) then
      ! In the current implementation, the Coulomb potential used for the non local potential is calculated in plane wave basis
      ! For details, please refer to Eq. 61 in doi:10.1016/j.cpc.2012.09.018
      call terminate_if_false( input%groundstate%Hybrid%BasisBareCoulomb == "pw", &
        "For RTTDDFT with hybrids only input%hybrid%barecoul%basis=pw is supported" )
      call set_barecoul_basis()
    end if
  end if

  if ( rt_inp%use_velocity_gauge() ) then
    if( rt_inp%pmat%read_pmat_from_file ) then 
      call terminate_if_false( file_pmat_exists( rt_inp%restart_file_handler, mpi_env_k ), &
        'File:'//trim( get_filename_pmat() )//' not found')
      call read_pmat( first_kpt, pmat, mpi_env_k, rt_inp%restart_file_handler )
      if ( molecular_dynamics%on ) then 
        call terminate_if_false( file_pmat_mt_exists( rt_inp%restart_file_handler, mpi_env_k ), &
          'File:'//trim( get_filename_pmat_mt() )//' not found')
        call read_pmat_mt( first_kpt, pmatmt, mpi_env_k, rt_inp%restart_file_handler )
      end if
    else
      if ( rt_inp%use_ks_basis() ) then
        call obtain_pmat_KSBasis( first_kpt, rt_inp%pmat%force_pmat_hermitian, apwalm, psi_gnd_lapwlo, pmat )
      else
        call obtain_pmat_LAPWloBasis( first_kpt, rt_inp%pmat%force_pmat_hermitian, apwalm, pmat, pmatmt )
      end if
    end if
    if( rt_inp%pmat%write_pmat_to_file ) then
      call write_pmat( first_kpt, pmat, mpi_env_k, rt_inp%restart_file_handler, kset_rttddft%nkpt )
      if ( molecular_dynamics%on ) call write_pmat_mt( first_kpt, pmatmt, mpi_env_k, rt_inp%restart_file_handler, kset_rttddft%nkpt )
    end if
  end if
  
  ! Attention: we allways need `IPA` to be false before the **first** call to `H%calculate`
  call H%set_IPA( .false. )
  call H%calculate( l_max_pot, apwalm, Gkset, psi_gnd_lapwlo )
  if ( rt_inp%use_ks_basis() .and. evolve_H0 ) then
    ! at t = 0, obtain the effective potential 
    call H%V_KS_0%copy_from( H%H_t ) 
    call H%V_KS_0%subtract( H%initial_eigenvalues )
  end if
  if ( .not. evolve_H0 .and. rt_inp%use_lapwlo_basis() ) call H%H_0%copy_from( H%H_t )
  ! Attention: now set `IPA` to its correct value
  call H%set_IPA( rt_inp%eeInteraction%use_ipa() )

  a_aux%components = 0._dp
  call overlap%initialize( apwalm, Gkset, pmatmt, a_aux, mathcalH=H%mathcalH, mathcalB=mathcalB )

  if ( rt_inp%use_velocity_gauge() ) then
    j_para_spurious%components = real_zero
    if ( rt_inp%subtract_J0 ) then
      call j_aux%evaluate_paramagnetic( psi, pmat, occupations, kset_rttddft%wkpt(first_kpt:last_kpt), mpi_env_k )
      j_para_spurious = j_aux%paramagnetic
    end if
  else
    e_aux%components = real_zero
    call get_td_overlap_det_and_berry_coupling_term( first_kpt, e_aux, pws_for_berry_phase, &
      psi, kset_rttddft, k_ptrs, td_overlap_det, berry_coupling_term )
    call p_vec_init%get_with_mtp( td_overlap_det, kset_rttddft%ngridk, kset_rttddft%ikmap, avec, prev_phases, .false. )
  end if

  if ( psi%has_frozen() ) then
    call update_density( first_kpt, psi, occupations, -1, .false., rt_inp%l_rad_step, &
      ks_lapwlo_transition_matrix=psi_gnd_lapwlo, dens_case=frozen )
    rhomt_frozen = rhomt
    rhoir_frozen = rhoir
  end if

  if( rt_inp%restart_previous_calculation() ) then
    call read_wavefunction( t, first_kpt, kset_rttddft%vkl(:, first_kpt:last_kpt), &
      psi%active, mpi_env_k, rt_inp%restart_file_handler )
    if( propagator%extrapolation_needed() ) then
      call read_wavefunction( t_minus_dt, first_kpt, kset_rttddft%vkl(:, first_kpt:last_kpt), &
        psi%active_save, mpi_env_k, rt_inp%restart_file_handler )
      call update_density( first_kpt, psi, occupations, 0, rt_inp%normalize_WF, &
        rt_inp%l_rad_step, rhomt_frozen, rhoir_frozen, psi_gnd_lapwlo, dens_case=save_and_frozen )
      call update_potential( coulomb_only =  rt_inp%eeInteraction%coulomb_only() )

      if ( rt_inp%use_lapwlo_basis() ) call overlap%calculate( apwalm, &
        Gkset, pmatmt, a_tot_t_minus_dt, mathcalH=H%mathcalH, mathcalB=mathcalB )

      if ( evolve_H0 ) call H%calculate( l_max_pot, apwalm, Gkset, psi_gnd_lapwlo, &
                              a_tot=a_tot_t_minus_dt, obtain_mathcalH=.true. )
      
      if ( rt_inp%use_velocity_gauge() ) then
        call H%add_external_coupling( a_tot_t_minus_dt, overlap, pmat )
      else
        call get_td_overlap_det_and_berry_coupling_term( first_kpt, e_vec_save, pws_for_berry_phase, psi, &
          kset_rttddft, k_ptrs, td_overlap_det, berry_coupling_term, .true. )
        call H%add_external_coupling( berry_coupling_term )
      end if
      call H%copy_H_t()
    end if ! propagator%extrapolation_needed()

    call update_density( first_kpt, psi, occupations, 0, rt_inp%normalize_WF, &
      rt_inp%l_rad_step, rhomt_frozen, rhoir_frozen, psi_gnd_lapwlo )
    call update_potential( coulomb_only = rt_inp%eeInteraction%coulomb_only() )

    if ( rt_inp%use_berry_phase() ) call read_phases( prev_phases, rt_inp%restart_file_handler, mpi_env_k )
  end if

  if( evolve_H0 .or. rt_inp%restart_previous_calculation() ) then
    call overlap%calculate( apwalm, Gkset, pmatmt, vec_pot%a_tot, mathcalH=H%mathcalH, mathcalB=mathcalB )
    call H%calculate( l_max_pot, apwalm, Gkset, psi_gnd_lapwlo, a_tot=vec_pot%a_tot, obtain_mathcalH=.true. )
    if ( rt_inp%use_velocity_gauge() ) then
      call H%add_external_coupling( vec_pot%a_tot, overlap, pmat )
    else
      call get_td_overlap_det_and_berry_coupling_term( first_kpt, e_vec, pws_for_berry_phase, &
        psi, kset_rttddft, k_ptrs, td_overlap_det, berry_coupling_term )
      call H%add_external_coupling( berry_coupling_term )
    end if
  end if

  if( rt_inp%do_from_scratch() .and. propagator%extrapolation_needed() ) call H%copy_H_t()

end subroutine

!> Allocate global MD arrays
subroutine allocate_MD_globals(first_kpt, last_kpt, allocate_mathcalB, allocate_B)
  !> Index of the first \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: first_kpt
  !> Index of the last \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: last_kpt
  !> if `.True`, we need to allocate the global array `mathcalB`
  logical, intent(in) :: allocate_mathcalB
  !> if `.True`, we need to allocate the global arrays `B_time` and `B_past`
  logical, intent(in) :: allocate_B

  if ( allocate_mathcalB ) allocate (mathcalB(nmatmax, nmatmax, 3, natmtot, first_kpt:last_kpt))
  if ( allocate_B ) then
    allocate ( B_time(nmatmax, nmatmax, first_kpt:last_kpt), source = zzero )
    allocate ( B_past(nmatmax, nmatmax, first_kpt:last_kpt), source = zzero )
  end if

end subroutine

!> Generate the basis type needed for the matrix elements evaluation and call `mt_init`
subroutine initialize_me( Gset )
  !> Set of \( \mathbf{G} \) vectors for LAPW expansion
  type(G_set), intent(in) :: Gset
  
  type(mt_basis_type) :: me_basis

  me_basis = mt_basis_type( spr(:, 1 : nspecies), nrmt(1 : nspecies), apwfr, lofr, &
    input%groundstate%lmaxapw, apword(:, 1 : nspecies), nlorb(1 : nspecies), lorbl(:, 1 : nspecies) )
  call me_init( me_basis, input%groundstate%lmaxvr, Gset )
  
end subroutine initialize_me

!> Output general information about the RT-TDDFT calculation using [[write_file_info]]
subroutine estimate_memory_and_write_to_info( ionDynamics, predictor_corrector, psi, &
    overlap, H, apwalm, kgrid_neighbours, n_kpt, pmat, pmatmt, &
    psi_gnd_lapwlo, pws_for_berry_phase, berry_coupling_term )
  !> Are we performing an MD calculation?
  logical, intent(in) :: ionDynamics
  !> if `.True`, the predictor corrector loop is employed
  logical, intent(in) :: predictor_corrector
  !> Basis-expansion coefficients of the KS-WFs
  class(wavefunction_set), intent(in) :: psi
  !> Overlap matrices
  class(overlap_set), intent(in) :: overlap
  !> Hamiltonian matrix
  class(hamiltonian_set), intent(in) :: H
  !> Matching coefficients of the (L)APWs
  complex(dp), intent(in) :: apwalm(:, :, :, :, :)
  !> Number of "neighbouring" \( \mathbf{k} \) points
  integer(i32), intent(in) :: kgrid_neighbours
  !> Total number of \( \mathbf{k} \) points
  integer(i32), intent(in) :: n_kpt
  !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
  !> (nmatmax, nmatmax, 3, first_kpt : last_kpt)
  complex(dp), intent(in), optional :: pmat(:, :, :, :)
  !> Muffin-tin part of the Momentum matrix
  !> (nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt)
  complex(dp), intent(in), optional :: pmatmt(:, :, :, :, :)
  !> KS-LAPW+lo transition matrix (ground state set in the LAPW+lo basis)
  complex(dp), intent(in), optional :: psi_gnd_lapwlo(:, :, :)
  !> Planewave matrix elements between neighbouring \( \mathbf{k} \) points
  complex(dp), intent(in), optional :: pws_for_berry_phase(:, :, :, :, :)
  !> Coupling with the external field through dynamical Berry phase approach
  complex(dp), intent(in), optional :: berry_coupling_term(:, :, :)

  character(len = str_128) :: string
  character(len=*), parameter :: formatMemory = '(A40,F12.1)'
  real(dp), parameter :: MB = 1048576._dp
  real(dp) :: aux_h, aux_w, aux_add, aux_td
  character(len=:), allocatable :: out_suffix

  aux_h = real( sizeof( overlap ) + sizeof( H ), dp ) / MB
  if( present( berry_coupling_term ) ) aux_h = aux_h + real( sizeof( berry_coupling_term ), dp ) / MB
  aux_w = real( sizeof( psi%frozen ) + sizeof( psi%active ) + sizeof( psi%groundstate ), dp ) / MB

  if ( kgrid_neighbours > 0 ) then
    aux_add = real( kgrid_neighbours, dp ) / real( kgrid_neighbours + psi%n_kpts(), dp ) &
    * real( sizeof( apwalm ) + sizeof( psi_gnd_lapwlo ), dp )
    ! Estimate memory to be allocated in [[get_td_overlap_det_and_berry_coupling_term]]
    aux_td = real( sizeof( psi%active ), dp ) * &
    real( n_kpt, dp ) / real( psi%n_kpts(), dp ) / MB + real( psi%n_occupied() * psi%n_occupied() * &
      psi%n_kpts() * 12._dp * sizeof( zzero ), dp ) / MB
  end if

  call write_file_info( 'Non-self-consistent GS for TDDFT calculations - finished' )
  call write_file_info_fill_line_with_char( '=' )
  call write_file_info( 'Allocated memory (MiB per MPI process)' )
  write (string, formatMemory) 'Coefficients to match LAPW functions:', real( sizeof( apwalm ), dp )/MB
  call write_file_info( string )
  write ( string, formatMemory ) 'Wavefunctions:', aux_w
  call write_file_info( string )
  write (string, formatMemory) 'Hamiltonian and Overlap matrices:', aux_h
  call write_file_info( string )
  if ( allocated( psi%active_save ) ) then
    write (string, formatMemory) 'Extra storage (WFs save):', &
      real( sizeof( psi%active_save ), dp ) / MB
    call write_file_info( string )
  end if
  if ( kgrid_neighbours > 0 ) then
    write (string, formatMemory) 'Extra storage (init):', aux_add / MB
    call write_file_info( string )
    write (string, formatMemory) 'To be allocated during propagation:', aux_td / MB
    call write_file_info( string )
  end if
  if ( predictor_corrector ) then
    write (string, formatMemory) 'Extra storage (predictor-corrector):', &
      real( sizeof( H%H_t%array ), dp ) / MB
    call write_file_info( string )
  end if
  if ( present ( psi_gnd_lapwlo ) ) then
    write (string, formatMemory) 'KS-LAPWlo transition matrix:', &
      real( sizeof( psi_gnd_lapwlo ), dp ) / MB
    call write_file_info( string )
  end if
  if ( present( pmat ) ) then
    write ( string, formatMemory ) 'Momentum matrix:', real( sizeof( pmat ), dp ) / MB
    call write_file_info( string )
  end if
  if ( ionDynamics ) then
    call write_file_info( string )
    write (string, formatMemory) 'Molecular Dynamics - Muffin-tin aux. matrices:', &
      real( sizeof(pmatmt) + sizeof(H%mathcalH) + sizeof(mathcalB) + sizeof(B_time) + sizeof(B_past), dp )/MB
    call write_file_info( string )
  end if
  if ( present( pws_for_berry_phase ) ) then
    write ( string, formatMemory ) 'PW MEs for dynamical Berry phase approach:', &
      real( sizeof( pws_for_berry_phase ), dp ) / MB
    call write_file_info( string )
  end if
  call write_file_info_fill_line_with_char('=')
  out_suffix = trim( filext )
  ! General info to be printed to RTTDDFT_INFO
  call write_file_info( 'Important output files: ' // filename_avec // out_suffix // ', ' &
                                                   // filename_pvec // out_suffix // ', ' &
                                                   // filename_jind // out_suffix // ', ' &
                                                   // filename_evec // out_suffix )
  call write_file_info( filename_jind // out_suffix // ' contains the x, y, and z components of the current density.' )
  call write_file_info( filename_pvec // out_suffix // ' contains the x, y, and z components of the polarization vector.' )
  call write_file_info( filename_avec // out_suffix // ' contains in each line 6 elements:' )
  call write_file_info( ': the x components of the induced and the total vector potential.' )
  call write_file_info( ': the y components of the induced and the total vector potential.' )
  call write_file_info( ': the z components of the induced and the total vector potential.' )
  call write_file_info( filename_evec // out_suffix // ' contains the x, y, and z components of the external electric field.' )
end subroutine

!> checks for consistency between gs hybrid calculation and rttddft and initializes pointer
subroutine adjustments_for_Hybrid_RTTDDFT()
  use rttddft_hybrids, only: hybrids_used
  logical :: is_compatible
  !> Tetrahedron method is used for Hybrid calculations
  call terminate_if_false(input%groundstate%stypenumber==-1, "stypenumber in the groundstate element must be set to libbzint to use To run RTTDDFT on top of hybrid functional calculations in the xs element")
  call terminate_if_false(hybrids_used(), "The non local current density can be computed only when hybrid functionals are used")
  is_compatible = .false.
  call is_gs_input_compatible_with_xs(input, is_compatible)
  if( is_compatible .eqv. .false.) call terminate_if_false( is_compatible, 'ERROR(rttddft_init): Parameters form GS not compatible with XS!')
  ! This pointer needs to get associated for when the calculation is started on top of a hybrid gs calculation
  if (.not. associated(input%groundstate%Hybrid)) &
      input%groundstate%Hybrid => getstructHybrid(emptynode)
end subroutine

!> For RT-TDDFT with hybrid funcionals, the same parameters for xs and the gs should be taken. This subroutine checks for that
subroutine is_gs_input_compatible_with_xs( inp, is_compatible)
  use modinput, only: input_type
  !> Information from the input file
  type(input_type), intent(in) :: inp
  !> True, if the input is compatible with rttddft and hybrid functionals (otherwise false)
  logical, intent(out) :: is_compatible
  !> Incompatibility message rttddft and hybrid functionals
  integer :: i

  is_compatible = .true.

  if (inp%xs%nosym .eqv. inp%groundstate%nosym) is_compatible = .false.
  do i = 1, 3
    if (inp%xs%ngridk(i) == inp%groundstate%ngridk(i)) is_compatible = .false.
    if (inp%xs%vkloff(i) == inp%groundstate%vkloff(i)) is_compatible = .false.
  end do
  if (inp%xs%reducek .eqv. inp%groundstate%reducek) is_compatible = .false.
  if (inp%xs%rgkmax == inp%groundstate%rgkmax) is_compatible = .false.
  if (inp%xs%swidth == inp%groundstate%swidth) is_compatible = .false.
  if (inp%xs%nempty == inp%groundstate%nempty) is_compatible = .false.

end subroutine

!> read WF and potential from potential gs run. For hybrid functionals, the parameters are read from the PBE run
subroutine read_WF_potential_rttddft( first_kpt, kset_rttddft, evecfv_gnd, occupations, initial_ks_energies )
  use modgw, only: kset
  !> Index of the first \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: first_kpt
  !> k set used in the RT-TDDFT module
  type(k_set), intent(in) :: kset_rttddft
  !> Basis-expansion coefficients of the groundstate KS-WFs
  complex(dp), intent(out) :: evecfv_gnd(:, :, first_kpt :)
  !> Initial occupations array
  real(dp), intent(out) :: occupations(:, first_kpt :)
  !> Initial eigenvalues array
  real(dp), intent(out) :: initial_ks_energies(:, first_kpt :)

  integer(i32) :: ik, last_kpt
  logical :: file_exists
  character(len=:), allocatable :: string

  last_kpt = ubound( evecfv_gnd, 3 )
  if ( hybrids_used() ) then
    inquire ( file='STATE_PBE.OUT', exist=file_exists )
    call terminate_if_false( file_exists, &
      'ERROR(rttddft_init): Start from GS calculation is not possible, STATE_PBE.OUT is missing!' )
    isreadstate0 = .false. ! We read not only from STATE.OUT
    string = filext
    filext = '_PBE.OUT'
  else
    string = filext
    filext = RTTDDFT_GND_sufix // filext
  end if
  call readstate()        ! read the density and potentials from file
  call gencore()          ! generate the core wavefunctions and densities
  call genmeffig()
  call linengy()          ! find the new linearization energies
  call genapwfr()         ! generate the APW radial functions
  call genlofr()          ! generate the local-orbital radial functions
  call olprad()           ! compute the overlap radial integrals
  if ( hybrids_used() ) then
    filext = string
    call energykncr()       ! core kinetic energy
    call init_product_basis()
    call readstate()
    call readfermi()
    call read_vxnl()        !The non local potential is read out from file
    call genmeffig()

    !----------------------------------------
    ! Read KS eigenvalues from file EVALFV.OUT
    !----------------------------------------
    if (allocated(evalfv)) deallocate (evalfv)
    allocate (evalfv(nstfv, kset%nkpt))
    evalfv(:, :) = 0.d0
    do ik = 1, kset%nkpt
      call getevalfv(kset%vkl(:, ik), evalfv(:, ik))
    end do

    ! VB / CB state index
    call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)

    ! The matrix sizes of mixed product basis quantities depend on if the core electrons are treated as valence
    ! This code block sets the dimension accordingly and is later read out in UpdateNonlocalCurrentDensity
    if ((input%gw%coreflag == 'all') .or. (input%gw%coreflag == 'xal')) then
      call Set_Dimension_mixed_product_basis(nomax + ncg)
    else
      call Set_Dimension_mixed_product_basis(nomax)
    end if

    ! Set BZ integration weights
    call kintw()

    deallocate (evalfv)

  end if

  ! Get the occupations, eigenvalues and eigenvectors from file
  do ik = first_kpt, last_kpt
    call getoccsv( kset_rttddft%vkl(:, ik), occupations(:, ik) )
  end do
  do ik = first_kpt, last_kpt
    call getevalsv( kset_rttddft%vkl(:, ik), initial_ks_energies(:, ik) )
  end do
  filext = string

  call read_wavefunction( groundstate, first_kpt, kset_rttddft%vkl(:, first_kpt:last_kpt), evecfv_gnd, mpi_env_k )

end subroutine


!> Get the initial KS states psi_gnd_lapwlo and basis set parameters apwalm
!> from the processes containing information about neighbouring \( \mathbf{k} \) points
subroutine get_data_from_neighbours( first_kpt, last_kpt, n_kpt, k_needed, &
  apwalm_extended, psi_gnd_lapwlo_extended, ik_to_array_position )
  !> Index of the first \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: first_kpt
  !> Last \( \mathbf{k} \) point for the current rank
  integer(i32), intent(in) :: last_kpt
  !> Total number of \( \mathbf{k} \) points
  integer(i32), intent(in) :: n_kpt
  !> Whether the current process needs \( \mathbf{k} \) point number i
  logical, intent(in) :: k_needed(:)
  !> apwalm containing info about host and neighbouring processes
  !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt + kgrid_neighbours)
  complex(dp), contiguous, intent(inout) :: apwalm_extended(:, :, :, :, first_kpt :)
  !> Initial states in LAPW basis from host and neighbouring processes
  !> (nmatmax, nstfv, first_kpt : last_kpt + kgrid_neighbours)
  complex(dp), contiguous, intent(inout) :: psi_gnd_lapwlo_extended(:, :, first_kpt :)
  !> Array mapping number of \( \mathbf{k} \) point (1 : nkpt) to the position in host array
  !> ( first_kpt : last_kpt + kgrid_neighbours )
  integer(i32), allocatable, intent(out) :: ik_to_array_position(:)

  integer(i32) :: ik, pos, broadcasting_rank, first_kpt_bcast, last_kpt_bcast
  complex(dp), allocatable :: apwlam_local(:, :, :, :), states_local(:, :)
  logical, allocatable :: k_needed_global(:)

  allocate ( ik_to_array_position( n_kpt ) )
  ik_to_array_position = -1
  do ik = first_kpt, last_kpt
    ik_to_array_position(ik) = ik 
  end do
  allocate( k_needed_global, source = k_needed )
  call xmpi_allreduce( k_needed_global, mpi_env_k )

  pos = last_kpt + 1
  allocate( apwlam_local, mold = apwalm_extended(:, :, :, :, first_kpt) )
  allocate( states_local, mold = psi_gnd_lapwlo_extended(:, :, first_kpt) )

  do broadcasting_rank = 0, mpi_env_k%procs - 1
    
    ! broadcasters' k points range
    first_kpt_bcast = firstk( broadcasting_rank, n_kpt )
    last_kpt_bcast = lastk( broadcasting_rank, n_kpt )

    do ik = first_kpt_bcast, last_kpt_bcast
      if ( .not. k_needed_global(ik) ) continue

      if ( rank == broadcasting_rank ) then
        apwlam_local = apwalm_extended(:, :, :, :, ik)
        states_local = psi_gnd_lapwlo_extended(:, :, ik)
      end if
      call xmpi_bcast( mpi_env_k, apwlam_local, broadcasting_rank )
      call xmpi_bcast( mpi_env_k, states_local, broadcasting_rank )    
   
      if ( k_needed(ik) ) then
        call assert( pos > ubound( apwalm_extended, 5 ), "more kgrid_neighbours found than exist" )

        apwalm_extended(:, :, :, :, pos) = apwlam_local
        psi_gnd_lapwlo_extended(:, :, pos) = states_local
        ik_to_array_position(ik) = pos
        pos = pos + 1
      end if
    end do
    call barrier( mpi_env_k )
  end do ! MPI procs

  pos = pos - 1
  call assert( pos /= ubound( apwalm_extended, 5 ), "less kgrid_neighbours found than exist" )
end subroutine get_data_from_neighbours


!> Evaluate the plane wave matrix elements for the neighbouring \( \mathbf{k} \) points
!> for each pair of the unperturbed KS states:
!> \[
!> W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{mn} = \int \psi_{\mathbf{k}}^* (\mathbf{r}) 
!> e^{-i \left(\mathbf{k}_{\alpha}^{\sigma} - \mathbf{k} \right) \mathbf{r}}  \psi_{\mathbf{k}_{\alpha}^{\sigma}} 
!> (\mathbf{r}) d \mathbf{r} = \int \psi_{\mathbf{k}}^* (\mathbf{r}) e^{-i \left(\mathbf{k}_{\alpha}^{\sigma} - \mathbf{k} \right) \mathbf{r}} 
!> \psi_{\mathbf{k}_{\alpha}^{\sigma} + s \mathbf{b}_{\alpha}} (\mathbf{r}) d \mathbf{r},
!> \]
!> where \( s = \pm 1 \), and \( \mathbf{b}_{\alpha} \) is the reciprocal lattice vector. 
subroutine calc_planewave_matrix_elements( first_kpt, k_ptrs, dk_vec, apwalm_extended, &
    psi_gnd_lapwlo_extended, pws_for_berry_phase, ik_to_array, Gkset, Gset, shift_positions, k_shifts )
  !> Index of the first \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: first_kpt
  !> Array containing indices of the neighbouring \( \mathbf{k} \) points
  integer(i32), intent(in) :: k_ptrs(:, :, :)
  !> Vectors dk in Cartesian coordinates for all directions and jumps
  real(dp), intent(in) :: dk_vec(:, :, :)
  !> Matching coefficients of the (L)APWs
  !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt + kgrid_neighbours)
  complex(dp), contiguous, intent(in) :: apwalm_extended(:, :, :, :, first_kpt :)
  !> (nmatmax, n_states, first_kpt : last_kpt + kgrid_neighbours)
  complex(dp), contiguous, intent(in) :: psi_gnd_lapwlo_extended(:, :, first_kpt :)
  !> Planewave matrix elements between neighbouring \( \mathbf{k} \) points
  !> (n_states, n_states, first_kpt : last_kpt, 3, 4)
  complex(dp), contiguous, intent(out) :: pws_for_berry_phase(:, :, first_kpt : , :, :)
  !> positions of local + neighbouring \( \mathbf{k} \) points in apwalm array
  integer(i32), intent(in) :: ik_to_array(:)
  !> set of \( \mathbf{G} + \mathbf{k} \) vectors for LAPW expansion
  type(Gk_set), intent(in) :: Gkset
  !> set of \( \mathbf{G} \) vectors for LAPW expansion
  type(G_set), intent(in) :: Gset
  !> Positions of the BZ shifts in the IR potential array
  integer(i32), intent(in) :: shift_positions(:)
  !> Positions of the shifts in the shift_positions array for each k point
  integer(i32), intent(in) :: k_shifts(:, :, :)

  integer(i32), parameter :: n_jumps = 4, n_cartesian_directions = 3, n_possible_shifts = 7
  integer(i32) :: n_states, k_left, k_right, k_direction, k_jump, last_kpt, is, ia, ias
  complex(dp), allocatable :: mt_contribution(:, :, :), plane_wave_sh(:, :), &
    mt_part(:, :, :, :, :)
  complex(dp) :: factor
  complex(dp), allocatable :: ir_pw(:, :), op_ir(:)

  ! get the arrrays' dimensions
  n_states = size( psi_gnd_lapwlo_extended, 2 )
  last_kpt = ubound( pws_for_berry_phase, 3 )

  pws_for_berry_phase = zzero

  call me_mt_alloc( mt_contribution )

  allocate( mt_part(size( mt_contribution, dim = 1 ), size( mt_contribution, dim = 2 ), &
    size( mt_contribution, dim = 3 ), n_cartesian_directions, n_jumps ), source = zzero )

  ! prepare MT PW expansion terms 4 \pi {\rm e}^{ -i \left(\bm{k}' - \bm{k} \right) Y_{l m} (\bm{k}) i^l j_l(kr)
  ! to sandwich between the MT basis functions
  do k_direction = 1, n_cartesian_directions
    do k_jump = 1, n_jumps

      if ( first_kpt == k_ptrs(first_kpt, k_direction, k_jump) ) continue

      mt_contribution = zzero
      do is = 1, nspecies
   
        call plane_wave_in_spherical_harmonics( -dk_vec(:, k_direction, k_jump), &
          spr(1 : nrmt(is), is), input%groundstate%lmaxvr, plane_wave_sh )
        
        do ia = 1, natoms(is)

          ias = idxas(ia, is)
          factor = exp( - zi * dot_product( dk_vec(:, k_direction, k_jump), atposc(:, ia, is) ) )

          ! compute gaunts times radial integrals
          call me_mt_prepare( is, ias, input%groundstate%lmaxvr, factor, &
            plane_wave_sh, zzero, mt_contribution(:, :, ias) )
            mt_part(:, :, ias, k_direction, k_jump) = mt_contribution(:, :, ias)

        end do ! ia
      end do ! is

    end do ! k_jump
  end do ! k_direction

  ! prepare the IR terms e^{ i s \bm{b}_{\alpha} \bm{r}} to sandwich between the PW basis functions
  call me_ir_alloc( ir_pw, n_possible_shifts )
  allocate( op_ir(Gset%ngrtot) )
  do is = 1, n_possible_shifts
    op_ir = zzero
    op_ir( shift_positions(is) ) = zone
    call zfftifc( n_cartesian_directions, Gset%ngrid, 1, op_ir )
    call me_ir_prepare( zone, op_ir, zzero, ir_pw(:, is) )
  end do

  do k_left = first_kpt, last_kpt
    do k_direction = 1, n_cartesian_directions
      do k_jump = 1, n_jumps

        k_right = k_ptrs(k_left, k_direction, k_jump)
        if ( k_right == k_left ) continue

        do is = 1, nspecies
          do ia = 1, natoms(is)

            ias = idxas(ia, is)

            call me_mt_mat( is, ias, ngk(1, k_left), ngk(1, k_right), &
              apwalm_extended(1 : ngk(1, k_left), :, :, ias, ik_to_array(k_left)), &
              apwalm_extended(1 : ngk(1, k_right), :, :, ias, ik_to_array(k_right)), &
              psi_gnd_lapwlo_extended(1 : nmat(1, k_left), :, ik_to_array(k_left)), &
              psi_gnd_lapwlo_extended(1 : nmat(1, k_right), :, ik_to_array(k_right)), &
              zone, mt_part(:, :, ias, k_direction, k_jump), zone, pws_for_berry_phase(:, :, k_left, k_direction, k_jump) )

          end do ! natoms
        end do ! nspecies

        ! ir contribution
        call me_ir_mat( Gkset, k_left, Gkset, k_right, &
          psi_gnd_lapwlo_extended(1 : nmat(1, k_left), :, ik_to_array(k_left)), &
          psi_gnd_lapwlo_extended(1 : nmat(1, k_right), :, ik_to_array(k_right)), &
          zone, ir_pw(:, k_shifts(k_left, k_direction, k_jump)), zone, &
          pws_for_berry_phase(:, :, k_left, k_direction, k_jump) )

      end do ! k_jump
    end do ! k_direction
  end do ! k_left

end subroutine calc_planewave_matrix_elements

!> Gather the information about neighbouring \( \mathbf{k} \) points based on the 
!> \( \mathbf{k} \) grid and MPI \( \mathbf{k} \)-sets
subroutine get_kgrid_neighbours_info( first_kpt, last_kpt, kset_rttddft, Gset, k_ptrs, dk_vec, &
    k_needed, proc_needed, kgrid_neighbours, shift_positions, k_shifts )
  !> Index of the first \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: first_kpt
  !> Index of the last \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: last_kpt
  !> Set of all \( \mathbf{k} \) points
  type(k_set), intent(in) :: kset_rttddft
  !> Set of \( \mathbf{G} \) vectors used for the matrix elements evaluation
  type(G_set), intent(in) :: Gset
  !> Array containing indices of the neighbouring \( \mathbf{k} \) points
  integer(i32), allocatable, intent(out) :: k_ptrs(:, :, :)
  !> Vectors dk in Cartesian coordinates for all directions and jumps
  real(dp), allocatable, intent(out) :: dk_vec(:, :, :)
  !> Array whose i's element is '.True.', if the current process needs information from \( \mathbf{k} \) point number i
  logical, allocatable, intent(out) :: k_needed(:)
  !> Array whose i's element is '.True.', if the current process needs information from the process number i - 1
  logical, allocatable, intent(out) :: proc_needed(:)
  !> Number of MPI processes containing \( \mathbf{k} \) points 'neighbouring' \( \mathbf{k} \) points from the current process
  integer(i32), intent(out) :: kgrid_neighbours
  !> Positions of the BZ shifts in the IR potential array
  integer(i32), allocatable, intent(out) :: shift_positions(:)
  !> Positions of the shifts in the shift_positions array for each k point
  integer(i32), allocatable, intent(out) :: k_shifts(:, :, :)

  integer(i32), parameter :: n_jumps = 4, n_cartesian_directions = 3, n_possible_shifts = 7
  integer(i32), parameter :: jumps(n_jumps) = [1, -1, 2, -2]
  integer(i32), parameter:: unit_vectors(n_cartesian_directions, n_cartesian_directions) = &
    reshape( [1, 0, 0, 0, 1, 0, 0, 0, 1 ], shape(unit_vectors) )
  integer(i32) :: ik, i, j, direction, iv(3), g(3)
  integer(i32), allocatable :: visited_k_pts(:), visited_procs(:), g_lattice(:, :)

  allocate( k_ptrs(kset_rttddft%nkpt, n_cartesian_directions, n_jumps) )
  k_ptrs = -1
  allocate( dk_vec(n_cartesian_directions, n_cartesian_directions, n_jumps ) )
  allocate( k_shifts(kset_rttddft%nkpt, n_cartesian_directions, n_jumps) )

  ! ik' = k_ptrs(ik, a, b), where a = 1, 2, 3 -- lattice direction
  ! b = +1, -1, +2, -2 -- \( \mathbf{k} \) grid step

  ! jumps to other BZ
  allocate( shift_positions(n_possible_shifts) )
  allocate( g_lattice(n_cartesian_directions, n_possible_shifts) )
  g_lattice = reshape( [0,0,0, &
                              1,0,0, &
                              -1,0,0, &
                              0,1,0, &
                              0,-1,0, &
                              0,0,1, &
                              0,0,-1], [n_cartesian_directions, n_possible_shifts] )
  do i = 1, n_possible_shifts
    g = g_lattice(:, i)
    shift_positions(i) = Gset%igfft( Gset%ivgig(g(1), g(2), g(3)) )
  end do

  do ik = 1, kset_rttddft%nkpt
    do direction = 1, n_cartesian_directions
      do j = 1, n_jumps
        iv = kset_rttddft%ivk(:, ik) + unit_vectors(:, direction) * jumps(j)
        if ( iv(direction) < 0 ) then
          k_shifts(ik, direction, j) = 2 * direction ! position in g_lattice
        else if ( iv(direction) >= kset_rttddft%ngridk(direction) ) then
          k_shifts(ik, direction, j) = 2 * direction + 1
        else
          k_shifts(ik, direction, j) = 1
        end if        
        iv(direction) = modulo( iv(direction), kset_rttddft%ngridk(direction) )
        k_ptrs(ik, direction, j) = kset_rttddft%ikmap(iv(1), iv(2), iv(3))
      end do
    end do
  end do
  do direction = 1, n_cartesian_directions
    do j = 1, n_jumps
      dk_vec(:, direction, j) = real( jumps(j), dp ) * kset_rttddft%bvec(:, direction ) &
        / real( kset_rttddft%ngridk(direction), dp )
    end do
  end do

  allocate( visited_k_pts(kset_rttddft%nkpt) )
  allocate( visited_procs(mpi_env_k%procs) )

  visited_k_pts = 0
  visited_procs = 0
  kgrid_neighbours = 0
  do i = 1, n_cartesian_directions
    do j = 1, n_jumps
      do ik = first_kpt, last_kpt
        if ( k_ptrs(ik, i, j) > last_kpt .or. k_ptrs(ik, i, j) < first_kpt ) then
          visited_k_pts(k_ptrs(ik, i, j)) = 1
          visited_procs(procofindex( k_ptrs(ik, i, j), kset_rttddft%nkpt, mpi_env_k%procs ) + 1) = 1
        end if
      end do
    end do
  end do
  
  kgrid_neighbours = sum( visited_k_pts )

  allocate( k_needed(kset_rttddft%nkpt), source = .false. )
  do i = 1, kset_rttddft%nkpt
    if ( visited_k_pts(i) == 1 ) k_needed(i) = .true.
  end do

  allocate( proc_needed(mpi_env_k%procs), source = .false. )
  do i = 1, mpi_env_k%procs
    if ( visited_procs(i) == 1 ) proc_needed(i) = .true.
  end do

end subroutine


!> Shift kpt_offset to the first parallelepiped if the input one is not there
!> for consistency with type(k_set)
subroutine ensure_valid_kpt_offset( kpt_offset )
  !> Shift of the \( \mathbf{k} \)-grid
  real(dp), intent(inout) :: kpt_offset(3)

  real(dp), parameter :: eps_lattice = 1.e-6_dp
  integer(i32) :: iv(3)

  if( any( abs( kpt_offset ) > 1.0_dp ) .or. any( kpt_offset < 0.0_dp ) ) then
    call r3frac( eps_lattice, kpt_offset, iv )
    call warning( 'Warning(initialize_rttddft): input%xs%vkloff mapped back to first k-parallelepiped: ' &
      // to_char( real( kpt_offset, sp ) ) ) 
  end if
end subroutine

!> Call the core exciting init routines
subroutine adjust_input_and_init_exciting_globals( global_input )
  !> Argument encapsulating input parameters
  type(input_type), intent(inout) :: global_input
  real(dp), parameter :: eps_rgkmax = 1.e-14_dp
  
  ! Backup groundstate variables
  call backup0()
  call backup1()
  if ( global_input%xs%rgkmax < eps_rgkmax ) &
    global_input%xs%rgkmax = global_input%groundstate%rgkmax
  call ensure_valid_kpt_offset( global_input%xs%vkloff )
  if ( hybrids_used() ) call adjustments_for_Hybrid_RTTDDFT()
  call mapxsparameters()
  call init0()
  call init1()
  call init2()
  if ( hybrids_used() ) call init_hybrids()

end subroutine

!> Get the energy gap using the KS energies and occupations array
subroutine get_energy_gap( first_kpt, ks_energies, occupations, n_kpt, energy_gap )
  !> Index of the first \( \mathbf{k} \) point treated by the current (MPI) rank
  integer(i32), intent(in) :: first_kpt
  !> Initial KS energies array (n_ks_states, n_kpt_this_proc)
  real(dp), contiguous, intent(in) :: ks_energies(:, first_kpt:)
  !> State occupations array (n_ks_states, n_kpt_this_proc)
  real(dp), contiguous, intent(in) :: occupations(:, first_kpt:)
  !> Total number of \( \mathbf{k} \) points
  integer(i32), intent(in) :: n_kpt
  !> Energy gap value 
  real(dp), intent(out) :: energy_gap

  real(dp), allocatable :: ks_energies_all_procs(:, :), occupations_all_procs(:, :)
  integer(i32) :: vbm_band_ind, cbm_band_ind, vbm_kpt_ind, cbm_kpt_ind, &
    gap_min_kpt_ind, last_kpt

  
  call assert( all( shape( ks_energies ) == shape( occupations ) ), "Incompatible energies &
    and occupations arrays provided to get_energy_gap." )
  last_kpt = ubound( ks_energies, 2 )

  allocate( ks_energies_all_procs( size( ks_energies, 1 ), n_kpt ), source = real_zero )
  ks_energies_all_procs(:, first_kpt : last_kpt) = ks_energies
  call xmpi_allgatherv( mpi_env_k, ks_energies_all_procs, size( ks_energies ) )

  allocate( occupations_all_procs, source = ks_energies_all_procs )
  occupations_all_procs(:, first_kpt : last_kpt) = occupations
  call xmpi_allgatherv( mpi_env_k, occupations_all_procs, size( occupations ) )

  call find_vbm_cbm( 1, size( ks_energies_all_procs, 1 ), size( ks_energies_all_procs, 2 ), &
    occupations_all_procs, ks_energies_all_procs, vbm_band_ind, cbm_band_ind, vbm_kpt_ind, &
    cbm_kpt_ind, gap_min_kpt_ind )

  energy_gap = ks_energies_all_procs( cbm_band_ind, cbm_kpt_ind ) - ks_energies_all_procs( vbm_band_ind, vbm_kpt_ind )
  if ( cbm_band_ind < vbm_band_ind ) energy_gap = 0._dp

end subroutine

end module rttddft_init
