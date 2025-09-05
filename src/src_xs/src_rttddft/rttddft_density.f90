! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created by Ronaldo Rodrigues Pela, May 2019
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module that manages what concerns charge density in RT-TDDFT calculations
module rttddft_Density
  use asserts, only: assert
  use precision, only: dp, i32
  use modmpi, only: mpi_env_k
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_density, timesec_RTTDDFT
  use constants, only: zzero, real_zero
  use mod_convergence, only : iscl
  use mod_potential_and_density, only: rhomt, rhoir
  use rttddft_Wavefunction, only: wavefunction_set
  use mod_rhoir, only: genrhoir
  use mod_rhovalk, only: rhovalk
  use general_matrix_multiplication, only: matrix_multiply
  use mod_eigensystem, only: nmat

  implicit none

  private

  public :: update_density, save_and_frozen, ground_state, frozen

  !> Enum with the density evaluation mode
  !> There are 4 options available: get density from `active` and `frozen` states, 
  !> `save` and `frozen` states, just `frozen` states, or all the unperturbed states.
  enum, bind(C)
    enumerator :: density_case
    enumerator :: active_and_frozen, save_and_frozen, frozen, ground_state
  end enum

contains
  !> In `update_density`, we obtain the charge density at time \(t\) 
  !> and put it in global `rho_mt` and `rho_ir` arrays. This routine 
  !> manages calls of `get_density_from_lapwlo_set` with appropriate arguments.
  subroutine update_density( first_kpt, psi, occupations, it, normalize, l_rad_step, &
      rhomt_frozen, rhoir_frozen, ks_lapwlo_transition_matrix, printTimings, t_dens, dens_case )
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Set of KS wavefunctions
    class(wavefunction_set), intent(in) :: psi
    !> Initial occupations array (n_states, n_kpts)
    real(dp), intent(in) :: occupations(:, :)
    !> Number of the current iteration (employed to give possible warnings)
    integer(i32), intent(in) :: it
    !> If `.true.`, normalize the charge density
    logical, intent(in) :: normalize
    !> Radial step length
    integer(i32), intent(in) :: l_rad_step
    !> Frozen part of the muffin-tin density (lmmaxvr, nrmtmax, natmtot)
    real(dp), optional, intent(in) :: rhomt_frozen(:, :, :)
    !> Frozen part of the IR density (ngrtot)
    real(dp), optional, intent(in) :: rhoir_frozen(:)
    !> KS-LAPW+lo transition matrix (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), optional, intent(in):: ks_lapwlo_transition_matrix(:, :, :)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(out) :: t_dens
    !> Enum telling which states should be used for density evaluation
    integer(kind( density_case )), optional, intent(in) :: dens_case

    integer(i32) :: last_kpt
    real(dp) :: ti 
    logical :: timings_general, timings_detailed, add_frozen
    integer(kind( density_case )) :: dens_case_
    real(dp), allocatable :: occupations_slice(:, :)
    complex(dp), allocatable :: local_lapw_set(:, :, :)

    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then
      call assert( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) call assert( timings_general, &
        'timings_general must be true if timings_detailed is true' )
    end if
    if( timings_general ) then
      call timesec( ti )
    end if

    dens_case_ = active_and_frozen
    if ( present( dens_case ) ) dens_case_ = dens_case
    add_frozen = .false. 
    if ( present( rhomt_frozen ) .or. present( rhoir_frozen ) ) then
      call assert( present( rhomt_frozen ) .and. present( rhoir_frozen ), &
        'Both contributions to frozen density should be provided to update_density' )
      ! for the ground state it is convenient to ignore precalculated frozen density
      add_frozen = ( .not. ( dens_case_ == ground_state ) )
    end if
    if ( dens_case_ == frozen ) then
      call assert( .not. add_frozen, &
        'frozen density requested with frozen rhoir and rhomt present in update_density' )
      call assert( psi%has_frozen(), &
        'frozen density requested with no frozen wavefunctions' )
    end if
    call assert( psi%n_kpts() == size( occupations, 2 ), 'psi and occupations have different n_kpts' )
    if ( .not. psi%expanded_in_lapwlo() ) call assert( present( ks_lapwlo_transition_matrix ), 'ks_lapwlo_transition_matrix is needed with KS basis' )

    rhomt = real_zero
    rhoir = real_zero
    last_kpt = first_kpt + psi%n_kpts() - 1
    select case ( dens_case_ )
    case( active_and_frozen )
      allocate( occupations_slice, source = occupations(psi%first_active() : psi%n_occupied(), :) )
      if ( .not. psi%expanded_in_lapwlo() ) call from_ks_to_lapw( ks_lapwlo_transition_matrix, nmat(1, first_kpt : last_kpt), &
        psi%active, local_lapw_set, printTimings, t_dens )
    case( save_and_frozen )
      allocate( occupations_slice, source = occupations(psi%first_active() : psi%n_occupied(), :) )
      if ( .not. psi%expanded_in_lapwlo() ) call from_ks_to_lapw( ks_lapwlo_transition_matrix, nmat(1, first_kpt : last_kpt), &
        psi%active_save, local_lapw_set, printTimings, t_dens )
    case( frozen )
      allocate( occupations_slice, source = occupations(1 : psi%n_frozen(), :) )
      if ( .not. psi%expanded_in_lapwlo() ) call from_ks_to_lapw( ks_lapwlo_transition_matrix, nmat(1, first_kpt : last_kpt), &
        psi%frozen, local_lapw_set, printTimings, t_dens )
    case( ground_state )
      allocate( occupations_slice, source = occupations )
      if ( .not. psi%expanded_in_lapwlo() ) call from_ks_to_lapw( ks_lapwlo_transition_matrix, nmat(1, first_kpt : last_kpt), &
        psi%groundstate, local_lapw_set, printTimings, t_dens )
    case default
      call assert( .false., 'unknown dens_case_' )
    end select

    if ( .not. psi%expanded_in_lapwlo() ) then
      call get_density_from_lapwlo_set( first_kpt, local_lapw_set, occupations_slice, it, &
        normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
    else
      select case ( dens_case_ )
      case( active_and_frozen )
        call get_density_from_lapwlo_set( first_kpt, psi%active, occupations_slice, it, &
          normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
      case( save_and_frozen )
        call get_density_from_lapwlo_set( first_kpt, psi%active_save, occupations_slice, it, &
          normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
      case( frozen )
        call get_density_from_lapwlo_set( first_kpt, psi%frozen, occupations_slice, it, &
          normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
      case( ground_state )
        call get_density_from_lapwlo_set( first_kpt, psi%groundstate, occupations_slice, it, &
          normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
      end select
    end if

    if( timings_general ) call timesec_RTTDDFT( ti, t_dens%total )

  end subroutine update_density

  !> Generate electron density from the set of KS states expanded over the LAPW+lo basis
  !> and put it into the `rho_mt` and `rho_ir` global arrays. The density is calculated
  !> using the same scheme as in the GS calcultations (refer to `scf_cycle.f90` 
  !> for the case `input%groundstate%useDensityMatrix`=`.false.`)
  subroutine get_density_from_lapwlo_set( first_kpt, psi_lapw_coeffs, occupations, it, &
    normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> States' expansion coefficients in the LAPW+lo set
    complex(dp), intent(in) :: psi_lapw_coeffs(:, :, :)
    !> Initial occupations array
    real(dp), intent(in) :: occupations(:, :)
    !> Number of the current iteration (employed to give possible warnings)
    integer(i32), intent(in) :: it
    !> If `.true.`, normalize the charge density
    logical, intent(in) :: normalize
    !> Radial step length
    integer(i32), intent(in) :: l_rad_step
    !> Whether the frozen part should be included
    logical, intent(in) :: add_frozen
    !> Frozen part of the muffin-tin density (lmmaxvr, nrmtmax, natmtot)
    real(dp), optional, intent(in) :: rhomt_frozen(:, :, :)
    !> Frozen part of the IR density (ngrtot)
    real(dp), optional, intent(in) :: rhoir_frozen(:)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(inout) :: t_dens

    logical :: timings_general, timings_detailed
    real(dp) :: ti
    integer(i32) :: i

    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then
      call assert( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) call assert( timings_general, &
        'timings_general must be true if timings_detailed is true' )
    end if
    if( timings_detailed ) call timesec( ti )

    call assert( size( psi_lapw_coeffs, 3 ) == size( occupations, 2 ), &
      "psi_lapw_coeffs and occupations have different n_kpts " )
    if ( add_frozen ) call assert( present( rhomt_frozen ) .and. present( rhoir_frozen ), &
      "rhomt_frozen and rhoir_frozen should be provided to get_density_from_lapwlo_set" )
    
    ! rhovalk and rhoir have omp critical inside
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(first_kpt, psi_lapw_coeffs, occupations) &
    !$OMP REDUCTION(+:rhomt, rhoir)
    do i = 1, size( occupations, 2 )
      call rhovalk( first_kpt + i - 1, psi_lapw_coeffs(:, :, i), occupations(:, i), rhomt )
      call genrhoir( first_kpt + i - 1, psi_lapw_coeffs(:, :, i), occupations(:, i), rhoir )
    end do
    !$OMP END PARALLEL DO

#ifdef MPI
    call mpisumrhoandmag( mpi_env_k )
#endif
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%rho )

    ! symmetrise the density
    call symrf( l_rad_step, rhomt, rhoir )
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%symrf )

    ! convert the density from a coarse to a fine radial mesh
    call rfmtctof( rhomt )
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%rfmtctof )

    ! generate the core wavefunctions and densities
    !call gencore

    if ( add_frozen ) then
      ! frozen density arrays contain core contribution along with the valence one
      rhoir = rhoir + rhoir_frozen
      rhomt = rhomt + rhomt_frozen
    else
      ! add the core density to the total density
      call addrhocr()
      if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%addrhocr )
    end if

    ! calculate the charges
    iscl = it
    call charge()
    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%charge )

    ! normalize the density
    if ( normalize ) then
      call rhonorm()
      if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%rhonorm )
    end if

  end subroutine

  !> Obtain the set of states' expansion coefficients \( {\bf C}_{j} \) in LAPW+lo basis
  !> from the expansion coefficients \( {\bf C}^{\rm KS}_{j} \)  in KS basis using
  !> the transition matrix \( M \):
  !> \[
  !> C_{ij} = \sum_{n} M_{in} C^{\rm KS}_{nj}.
  !> \]
  subroutine from_ks_to_lapw( transition_matrix, k_dependent_basis_size, set_in_ks_basis, &
      set_in_lapw_basis, printTimings, t_dens )
    !> Transition matrix (change-of-basis matrix) \( M \) used to change the basis from 
    !> the KS one to the LAPW+lo set (n_basis_lapw_max, n_basis_ks, n_kpts)
    complex(dp), intent(in) :: transition_matrix(:, :, :)
    !> k-dependent lapw basis dimensions array (n_kpts)
    integer(i32), intent(in) :: k_dependent_basis_size(:)
    !> KS basis expansion coefficients \( {\bf C}^{\rm KS}_{j} \) (n_basis_ks, n_states, n_kpts)
    complex(dp), intent(in) :: set_in_ks_basis(:, :, :)
    !> LAPW+lo basis expansion coefficients \( {\bf C}_{j} \) (n_basis_lapw_max, n_states, n_kpts)
    complex(dp), allocatable, intent(out) :: set_in_lapw_basis(:, :, :)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(inout) :: t_dens

    integer(i32) :: n_states, ik, n_kpts
    logical :: timings_general, timings_detailed
    real(dp) :: ti

    n_states = size( set_in_ks_basis, 2 )
    n_kpts = size( transition_matrix, 3 )

    call assert( n_kpts == size( k_dependent_basis_size ), &
      'transition_matrix and k_dependent_basis_size have different n_kpts' )
    call assert( n_kpts == size( set_in_ks_basis, 3 ), &
      'transition_matrix and set_in_ks_basis have different n_kpts' )
    call assert( size( transition_matrix, 2 ) == size( set_in_ks_basis, 1 ), &
      'transition_matrix and set_in_ks_basis have different n_basis' )
    
    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then
      call assert( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) call assert( timings_general, &
        'timings_general must be true if timings_detailed is true' )
    end if
    if( timings_detailed ) call timesec( ti )

    allocate( set_in_lapw_basis( size( transition_matrix, 1 ), n_states, n_kpts ), source = zzero )

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik), &
    !$OMP SHARED(n_kpts, transition_matrix, k_dependent_basis_size, n_states) &
    !$OMP SHARED(set_in_ks_basis, set_in_lapw_basis)
    !$OMP DO
    do ik = 1, n_kpts
      call matrix_multiply( transition_matrix( 1 : k_dependent_basis_size(ik), :, ik), &
        set_in_ks_basis(:, :, ik), set_in_lapw_basis(1 : k_dependent_basis_size(ik), :, ik) )
    end do
    !$OMP END PARALLEL

    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%basis )

  end subroutine

end module rttddft_Density
