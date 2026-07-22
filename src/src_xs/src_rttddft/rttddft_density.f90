!> Module that manages what concerns charge density in RT-TDDFT calculations
module rttddft_Density
#include "asserts.fpp"
  use constants, only: zzero, real_zero
  use general_matrix_multiplication, only: matrix_multiply
  use mod_potential_and_density, only: magir, magmt, rhoir, rhomt
  use mod_rhoir, only: genrhoir
  use mod_rhovalk, only: rhovalk
  use modmpi, only: mpi_env_k
  use precision, only: dp, i32
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_density, timesec_RTTDDFT
  use rttddft_Wavefunction, only: wavefunction_set
  use to_char_conversion, only: to_char

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
  subroutine update_density( psi, it, normalize, l_rad_step, &
      rhomt_frozen, rhoir_frozen, printTimings, t_dens, dens_case )
    !> Set of KS wavefunctions
    class(wavefunction_set), intent(in) :: psi
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
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the electronic density
    type(Timing_RTTDDFT_density), optional, intent(out) :: t_dens
    !> Enum telling which states should be used for density evaluation
    integer(kind( density_case )), optional, intent(in) :: dens_case

    integer(i32) :: first_kpt, last_kpt, n1, n2, ik
    real(dp) :: ti 
    logical :: timings_general, timings_detailed, add_frozen
    integer(kind( density_case )) :: dens_case_
    real(dp), allocatable :: occupations_slice(:, :)
    complex(dp), allocatable :: local_lapw_set(:, :, :), second_variation_set(:, :, :)

    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then
      CALL_ASSERT( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) then
        CALL_ASSERT( timings_general,  'timings_general must be true if timings_detailed is true' )
      end if
    end if
    if( timings_general ) then
      call timesec( ti )
    end if

    dens_case_ = active_and_frozen
    if ( present( dens_case ) ) dens_case_ = dens_case
    add_frozen = .false. 
    if ( present( rhomt_frozen ) .or. present( rhoir_frozen ) ) then
      CALL_ASSERT( present( rhomt_frozen ) .and. present( rhoir_frozen ),  'Both contributions to frozen density should be provided to update_density' )
      ! for ground state, ignore precalculated frozen density
      add_frozen = ( .not. ( dens_case_ == ground_state ) )
    end if
    if ( dens_case_ == frozen ) then
      CALL_ASSERT( .not. add_frozen,  'frozen density requested with frozen rhoir and rhomt present in update_density' )
      CALL_ASSERT( psi%has_frozen(),  'frozen density requested with no frozen wavefunctions' )
    end if

    rhomt = real_zero
    rhoir = real_zero
    if( psi%has_second_variation() ) then
      magir = real_zero
      magmt = real_zero
    end if
    first_kpt = psi%first_kpt()
    last_kpt = psi%last_kpt()
    select case ( dens_case_ )
    case( active_and_frozen )
      n1 = psi%first_active(); n2 = psi%n_occupied()
      if ( psi%expanded_in_ks() ) then
        if( psi%has_second_variation() ) then
          local_lapw_set = psi%groundstate_lapwlo
          allocate( second_variation_set, mold = psi%active )
          do ik = first_kpt, last_kpt
            call matrix_multiply( psi%groundstate_second_variation(:, :, ik), &
              psi%active(:, :, ik), second_variation_set(:, :, ik) )
          end do
        else
          call wrapper_from_ks_to_lapw( psi%active )
        end if
      end if
    case( save_and_frozen )
      CALL_ASSERT( .not. psi%has_second_variation(), 'second variation not supported when frozen states are present' )
      n1 = psi%first_active(); n2 = psi%n_occupied()
      if ( psi%expanded_in_ks() ) call wrapper_from_ks_to_lapw( psi%active_save )
    case( frozen )
      CALL_ASSERT( .not. psi%has_second_variation(), 'second variation not supported when frozen states are present' )
      n1 = 1; n2 = psi%n_frozen()
      if ( psi%expanded_in_ks() ) call wrapper_from_ks_to_lapw( psi%frozen )
    case( ground_state )
      n1 = 1; n2 = psi%n_occupied()
      if ( psi%expanded_in_ks() ) local_lapw_set = psi%groundstate_lapwlo
      if ( psi%has_second_variation() ) second_variation_set = psi%groundstate_second_variation
    case default
      CALL_ASSERT( .false., 'unknown dens_case_' )
    end select
    occupations_slice = psi%occupations(n1:n2, :)

    if ( psi%expanded_in_ks() ) then
      call wrapper_get_density_from_lapwlo_set( local_lapw_set )
    else
      select case ( dens_case_ )
      case( active_and_frozen )
        call wrapper_get_density_from_lapwlo_set( psi%active )
      case( save_and_frozen )
        call wrapper_get_density_from_lapwlo_set( psi%active_save )
      case( frozen )
        call wrapper_get_density_from_lapwlo_set( psi%frozen )
      case( ground_state )
        call wrapper_get_density_from_lapwlo_set( psi%groundstate_lapwlo )
      end select
    end if

    if( timings_general ) call timesec_RTTDDFT( ti, t_dens%total )

    contains
    subroutine wrapper_get_density_from_lapwlo_set( psi_lapw )
      complex(dp), contiguous, intent(in) :: psi_lapw(:, :, :)

      ! The following code should be ok. However gfortran (v. 14) gives a seg. fault
      ! call get_density_from_lapwlo_set( first_kpt, psi_lapw, second_variation_set, &
      !     occupations_slice, it, normalize, l_rad_step, add_frozen, rhomt_frozen=rhomt_frozen, rhoir_frozen=rhoir_frozen, &
      !     printTimings=printTimings, t_dens=t_dens )
      if( present( rhomt_frozen ) ) then
        call get_density_from_lapwlo_set( first_kpt, psi_lapw, second_variation_set, &
          occupations_slice, it, normalize, l_rad_step, add_frozen, rhomt_frozen, &
          rhoir_frozen, printTimings, t_dens )
      else
        call get_density_from_lapwlo_set( first_kpt, psi_lapw, second_variation_set, &
          occupations_slice, it, normalize, l_rad_step, add_frozen, printTimings=printTimings, t_dens=t_dens )
      end if

    end subroutine

    subroutine wrapper_from_ks_to_lapw( psi_ks )
      complex(dp), contiguous, intent(in) :: psi_ks(:, :, :)

      call from_ks_to_lapw( psi%groundstate_lapwlo, psi_ks, local_lapw_set, printTimings, t_dens )
    end subroutine

  end subroutine update_density

  !> Generate electron density from the set of KS states expanded over the LAPW+lo basis
  !> and put it into the `rho_mt` and `rho_ir` global arrays. The density is calculated
  !> using the same scheme as in the GS calcultations (refer to `scf_cycle.f90` 
  !> for the case `input%groundstate%useDensityMatrix`=`.false.`)
  subroutine get_density_from_lapwlo_set( first_kpt, psi_lapw_coeffs, psi_second_variation, occupations, it, &
    normalize, l_rad_step, add_frozen, rhomt_frozen, rhoir_frozen, printTimings, t_dens )
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> States' expansion coefficients in the LAPW+lo set
    complex(dp), contiguous, intent(in) :: psi_lapw_coeffs(:, :, :)
    !> Second variation coefficients
    complex(dp), contiguous, optional, intent(in) :: psi_second_variation(:, :, :)
    !> Initial occupations array
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Number of the current iteration (employed to give possible warnings)
    integer(i32), intent(in) :: it
    !> If `.true.`, normalize the charge density
    logical, intent(in) :: normalize
    !> Radial step length
    integer(i32), intent(in) :: l_rad_step
    !> Whether the frozen part should be included
    logical, intent(in) :: add_frozen
    !> Frozen part of the muffin-tin density (lmmaxvr, nrmtmax, natmtot)
    real(dp), contiguous, optional, intent(in) :: rhomt_frozen(:, :, :)
    !> Frozen part of the IR density (ngrtot)
    real(dp), contiguous, optional, intent(in) :: rhoir_frozen(:)
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
      CALL_ASSERT( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) then
        CALL_ASSERT( timings_general, 'timings_general must be true if timings_detailed is true' )
      end if
    end if
    if( timings_detailed ) call timesec( ti )

    CALL_ASSERT( size( psi_lapw_coeffs, 3 ) == size( occupations, 2 ), "psi_lapw_coeffs and occupations have different n_kpts " )
    if ( add_frozen ) then 
      CALL_ASSERT( present( rhomt_frozen ) .and. present( rhoir_frozen ), "rhomt_frozen and rhoir_frozen absent" )
    end if
    
    if( present( psi_second_variation ) ) then
      CALL_ASSERT( size( psi_second_variation, 3 ) == size( occupations, 2 ), "psi_second_variation and occupations have different n_kpts " )
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(first_kpt, psi_lapw_coeffs, occupations, psi_second_variation) &
      !$OMP REDUCTION(+:rhomt, rhoir, magmt, magir)
      do i = 1, size( occupations, 2 )
        call rhovalk( first_kpt + i - 1, psi_lapw_coeffs(:, :, i), occupations(:, i), &
          rhomt, magmt, psi_second_variation(:, :, i) )
        call genrhoir( first_kpt + i - 1, psi_lapw_coeffs(:, :, i), occupations(:, i), &
          rhoir, magir, psi_second_variation(:, :, i) )
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(first_kpt, psi_lapw_coeffs, occupations) &
      !$OMP REDUCTION(+:rhomt, rhoir)
      do i = 1, size( occupations, 2 )
        call rhovalk( first_kpt + i - 1, psi_lapw_coeffs(:, :, i), occupations(:, i), rhomt )
        call genrhoir( first_kpt + i - 1, psi_lapw_coeffs(:, :, i), occupations(:, i), rhoir )
      end do
      !$OMP END PARALLEL DO
    end if

    call mpisumrhoandmag( mpi_env_k )
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
    call charge( 'real-time propagation step ' // to_char( it ) )
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
  subroutine from_ks_to_lapw( transition_matrix, set_in_ks_basis, &
      set_in_lapw_basis, printTimings, t_dens )
    !> Transition matrix (change-of-basis matrix) \( M \) used to change the basis from 
    !> the KS one to the LAPW+lo set (n_basis_lapw_max, n_basis_ks, n_kpts)
    complex(dp), intent(in) :: transition_matrix(:, :, :)
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

    CALL_ASSERT( n_kpts == size( set_in_ks_basis, 3 ), 'transition_matrix and set_in_ks_basis have different n_kpts' )
    CALL_ASSERT( size( transition_matrix, 2 ) == size( set_in_ks_basis, 1 ), 'transition_matrix and set_in_ks_basis have different n_basis' )
    
    timings_general = .false.
    timings_detailed = .false.
    if( present( printTimings ) ) then
      CALL_ASSERT( present( t_dens ), 'Optional argument t_dens must also be present' )
      call printTimings%get( timings_general, timings_detailed )
      if( timings_detailed ) then
        CALL_ASSERT( timings_general,  'timings_general must be true if timings_detailed is true' )
      end if
    end if
    if( timings_detailed ) call timesec( ti )

    allocate( set_in_lapw_basis( size( transition_matrix, 1 ), n_states, n_kpts ), source = zzero )

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik), &
    !$OMP SHARED(n_kpts, transition_matrix, n_states) &
    !$OMP SHARED(set_in_ks_basis, set_in_lapw_basis)
    !$OMP DO
    do ik = 1, n_kpts
      call matrix_multiply( transition_matrix(:, :, ik), &
        set_in_ks_basis(:, :, ik), set_in_lapw_basis(:, :, ik) )
    end do
    !$OMP END PARALLEL

    if( timings_detailed ) call timesec_RTTDDFT( ti, t_dens%basis )

  end subroutine

end module rttddft_Density
