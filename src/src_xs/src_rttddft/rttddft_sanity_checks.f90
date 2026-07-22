!> This module handles checks of the user-provided input for the RT module. The calculations are 
!> terminated, if the input is inconsistent or incompatible with the current version of the code.
module rttddft_sanity_checks
#include "asserts.fpp"
  use constants, only: real_zero
  use exciting_mpi, only: xmpi_allgatherv
  use general_find_vbm_cbm, only: find_vbm_cbm
  use math_utils, only: all_zero
  use mod_mpi_env, only: mpiinfo
  use modinput, only: input_type, isspinspiral
  use modmpi, only: terminate, terminate_if_false
  use physical_constants, only: c
  use precision, only: dp, i32, sp
  use rttddft_electric_field, only: Electric_Field
  use rttddft_VectorPotential, only: Vector_Potential
  use rttddft_Wavefunction, only: wavefunction_set
  use rttddft_Wavefunction, only: wavefunction_set
  use to_char_conversion, only: to_char
  use vector_multiplication, only: norm

  implicit none
  
  private
  public :: check_rttddft_input, check_rttddft_setup

  integer(i32), parameter :: n_cartesian_directions = 3
  real(dp), parameter :: &
    e_field_squared_au_to_intensity_wcm2 = 3.50941e16_dp, &
    intensity_extremely_high = 1.e18_dp, &
    eps_kick_width = 1.e-14_dp, &
    t_step_scale = 0.2_dp, &
    eps_energy_gap = 1.e-5_dp, &
    eps_e_field = 1.e-12_dp, &
    eps_scissor = 1.e-10_dp, &
    e_field_extremely_high = sqrt( intensity_extremely_high / e_field_squared_au_to_intensity_wcm2 )
  
contains

  !> Check if variables given in the input file make sense
  subroutine check_rttddft_input( inp )
    !> type with the variables given in the input file
    type(input_type), intent(in) :: inp

    integer(i32) :: i
    
    call terminate_if_false( .not. inp%groundstate%solver%packedmatrixstorage, &
      & 'RT-TDDFT does not work with matrices stored in a packed form.' )

    ! Consistency check for spin-polarized calculations
    if( inp%groundstate%tevecsv .or. associated( inp%groundstate%spin ) ) then
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%fieldCoupling ) == "velocityGauge", &
        "Only velocity gauge is currently available for spin-polarised calculations." )
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%basis ) == "unperturbedKS", &
        "Only unperturbedKS is currently available for spin-polarised calculations." )
      if( isspinspiral() ) call terminate( "RT-TDDFT not implemented for spin-spiral calculations." )
      call terminate_if_false( inp%xs%realTimeTDDFT%numberOfFrozenStates == 0, &
        "No state freezing currently possible for spin-polarised calculations." )
      call terminate_if_false( inp%groundstate%spin%realspace, &
        "RT-TDDFT only implemented for realspace spin-polarised calculations." )
      call terminate_if_false( all_zero( inp%groundstate%spin%bfieldc ), &
        "RT-TDDFT only implemented when bfieldc is zero." )
      call terminate_if_false( associated(inp%xs%realTimeTDDFT%spinPropagation), &
        "spin Propagation must be present for spin-polarized calculations." )
      ! all checks for SOC
      if( inp%xs%realTimeTDDFT%spinPropagation%SOCGaugeCorrections ) &
        call terminate( "SOC gauge corrections not yet implemented." )
      if( associated( inp%groundstate%spin ) ) then
        if( inp%groundstate%spin%spinorb ) then 
          if( .not. inp%xs%realTimeTDDFT%spinPropagation%updateSOC ) then
            call warn( "updateSOC=false may lead to inaccurate results." )
          else
            if( .not. inp%xs%realTimeTDDFT%spinPropagation%SOCGaugeCorrections ) &
              call warn( "SOCgaugeCorrections=false may lead to unphysical results, as it breaks gauge invariance." )
          end if
        end if
      end if
    end if
    ! iora*
    if( trim(inp%groundstate%ValenceRelativity) == "iora*" ) call terminate( &
      & 'RT-TDDFT not implemented for ValenceRelativity="iora*"'   )

    ! Consistency check: laser has been defined?
    call terminate_if_false( associated( inp%xs%realTimeTDDFT%laser ), &
      & 'Element <laser> in <realTimeTDDFT> not found')

    if( associated(inp%xs%realTimeTDDFT%predictorCorrector) ) then
      ! Consistency check: MD and predictor corrector?
      call terminate_if_false( .not. associated( inp%MD ), &
        & 'It is currently not possible to use the predictor corrector method together with molecular dynamics' )
      ! Consistency check: predictor corrector method cannot be used with propagators SE and EH
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%propagator ) /= 'SE' .and. &
        trim( inp%xs%realTimeTDDFT%propagator ) /= 'EH', 'EH and SE methods are not compatible with predictor-corrector' )
      ! Consistency check: predictor corrector method should not be used with frozen ee interaction
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%eeInteraction ) /= "IPA", &
        & 'Predictor corrector method should not be used together with IP approximation')
    end if

    if ( trim( inp%xs%realTimeTDDFT%fieldCoupling ) == "berryPhase" ) then
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%basis ) == "unperturbedKS", &
        "Berry-phase coupling is currently available only with the KS basis" )
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%laser%fieldType ) == "total", &
        "Berry-phase coupling is currently available only with total field given" )
    else
      call terminate_if_false( associated( inp%xs%realTimeTDDFT%pmat ), &
      & 'Element <pmat> in <realTimeTDDFT> not found' )
    end if

    if( associated( inp%MD ) ) then
      if( trim( inp%xs%realTimeTDDFT%do ) /= "fromscratch" ) then 
        call terminate_if_false( trim( inp%xs%realTimeTDDFT%propagator ) == "SE" .or. trim( inp%xs%realTimeTDDFT%propagator ) == "EH", &
          "Restart for Ehrenfest MD is currently only implemented for the SE and EH propagators" )
        call terminate_if_false( inp%xs%realTimeTDDFT%timeStep == inp%MD%timeStep, &
          "Restart for Ehrenfest MD is currently only implemented when the RT-TDDFT and MD timesteps are the same")
      end if
      call terminate_if_false( inp%xs%realTimeTDDFT%numberOfFrozenStates == 0, &
        "No state freezing currently possible for MD calculations" )
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%basis ) == "LAPWlo", &
        "Usage of the KS basis set is currently unavailable for MD calculations" )
      call terminate_if_false( abs( inp%xs%scissor ) < eps_scissor, &
        "Scissor correction is currently unavailable for MD calculations" )
    end if

    if ( inp%xs%realTimeTDDFT%calculateTotalEnergy ) then
      ! Consistency check: real-time total energy is ill-defined with frozen ee interaction
      call terminate_if_false( trim( inp%xs%realTimeTDDFT%eeInteraction ) /= "IPA", &
        & 'Real-time total energy should not be evaluated with IP approximation')
    end if

    if ( associated( inp%xs%realTimeTDDFT%laser%kickarray ) ) then
      associate( kick_array => inp%xs%realTimeTDDFT%laser%kickarray )
        do i = 1, size( kick_array )
          if ( abs( kick_array(i)%kick%width ) < eps_kick_width ) &
            call warn( 'electric field is ill-defined at time ' // &
              to_char( real( kick_array(i)%kick%t0, sp) ) // ' for the kick number ' // to_char( i ) )
        end do
      end associate
    end if

    if ( inp%xs%realTimeTDDFT%numberOfFrozenStates > 0 .and. trim( inp%xs%realTimeTDDFT%fieldCoupling ) == "velocityGauge"  ) then
      if ( .not. inp%xs%realTimeTDDFT%orthogonalizeAgainstFrozen ) call warn( &
        'It is strongly recommended to set orthogonalizeAgainstFrozen="true" with velocity gauge if frozen states are present.' )
      if ( .not. inp%xs%realTimeTDDFT%subtractJ0 ) call warn( &
        'It is recommended to set subtractJ0="true" with velocity gauge if frozen states are present.'  )
    end if

    contains
      subroutine warn( message )
        character(len=*), intent(in) :: message
        character(len=*), parameter :: procedure_name = "check_rttddft_input"
        character(len=*), parameter :: warning_header = "Warning(" // procedure_name // "): "

        call wrapper( warning_header, message )
      end subroutine

  end subroutine

  ! Check if input variables make sense after the RT module initialization
  subroutine check_rttddft_setup( mpi_env, time_step, initial_ks_energies, use_berry_phase, &
    vec_pot, t_start, t_end, lattice_vectors, psi, scissor )
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Time evolution step
    real(dp), intent(in) :: time_step
    !> Initial KS energies array (n_ks_states, n_kpt_current_rank)
    real(dp), contiguous, intent(in) :: initial_ks_energies(:, :)
    !> Whether the Berry-phase coupling is used
    logical, intent(in) :: use_berry_phase
    !> Argument that encapsulates the vector potential
    type(Vector_Potential), intent(in) :: vec_pot
    !> Time of the start of the evolution
    real(dp), intent(in) :: t_start
    !> Time of the end of the evolution
    real(dp), intent(in) :: t_end
    !> Array containing the lattice vectors
    real(dp), intent(in) :: lattice_vectors(:, :)
    !> Set of KS states
    class(wavefunction_set), intent(in) :: psi
    !> Requested value of a scissor correction
    real(dp), intent(in) :: scissor

    real(dp) :: energy_gap, t_step_critical, e_field_critical(n_cartesian_directions), &
      e_field_max_lattice(n_cartesian_directions), lattice_vector_norm(n_cartesian_directions), &
      e_field_max_magnitude
    type(Electric_Field) :: e_aux
    integer(i32) :: i, j

    call get_energy_gap( mpi_env, initial_ks_energies, psi%occupations, energy_gap )
    if ( scissor > eps_scissor ) call terminate_if_false( energy_gap > eps_energy_gap, &
      "Scissor correction requested for a non-insulating material." )

    do i = 1, n_cartesian_directions
      lattice_vector_norm(i) = norm( lattice_vectors(:, i) )
    end do

    e_field_max_lattice = real_zero
    e_field_max_magnitude = real_zero
    do j = 1, int( (t_end - t_start) / time_step, kind = i32 )
      e_aux%components = - vec_pot%get_dA_dt( t_start + real( j, dp ) * time_step ) / c
      do i = 1, n_cartesian_directions
        e_field_max_lattice(i) = max( e_field_max_lattice(i), &
          abs( dot_product( e_aux%components, lattice_vectors(:, i) ) ) / lattice_vector_norm(i) )
        e_field_max_magnitude = max( e_field_max_magnitude, norm( e_aux%components ) )
      end do
    end do

    if ( all( e_field_max_lattice < eps_e_field ) ) call warn( "external field amplitude is zero." )

    if ( use_berry_phase ) then
      call terminate_if_false( energy_gap > eps_energy_gap, " &
        Berry-phase field coupling is only defined for an insulator." )
      do i = 1, n_cartesian_directions
        if ( e_field_max_lattice(i) > eps_e_field ) call terminate_if_false( psi%kset%ngridk(i) > 2, &
          "At least 3 k-points in lattice direction " // to_char(i) // " are needed for the &
          Berry-phase coupling operator construction." )

        e_field_critical(i) = energy_gap / ( real( psi%kset%ngridk(i), dp ) * &
          sqrt( dot_product( lattice_vectors(:, i), lattice_vectors(:, i) ) ) )
        if ( e_field_max_lattice(i) > e_field_critical(i) ) call warn( &
          'field strength ' // to_char( real( e_field_max_lattice(i), sp) ) // ' in lattice direction ' &
          // to_char(i) // ' exceeds the estimated largest reasonable value of ' &
          // to_char( real( e_field_critical(i), sp) ) // ', see Zener &
          tunneling discussion in [PRL 89, 117602 (2002), PRB 69, 085106 (2004)].' )
      end do
    end if
    
    if ( e_field_max_magnitude > e_field_extremely_high ) call warn( &
      'field strength magnitude corresponds to extremely high laser intensity of ' &
      // to_char( real( e_field_squared_au_to_intensity_wcm2 * e_field_max_magnitude**2, sp) ) // ' W/cm^2.' )

    t_step_critical = t_step_scale / &
    ( maxval( initial_ks_energies ) - minval( initial_ks_energies ) )
    if ( time_step > t_step_critical ) call warn( &
      'time step ' // to_char( real( time_step, sp ) ) // ' exceeds the roughly-estimated &
      largest reasonable value of ' // to_char( real( t_step_critical, sp) ) // '.' )

    contains
      subroutine warn( message )
        character(len=*), intent(in) :: message
        character(len=*), parameter :: procedure_name = "check_rttddft_setup"
        character(len=*), parameter :: warning_header = "Warning(" // procedure_name // "): "

        call wrapper( warning_header, message )
      end subroutine
  end subroutine

  !> (private) Get the energy gap using the KS energies and occupations array
  subroutine get_energy_gap( mpi_env, ks_energies, occupations, energy_gap )
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Initial KS energies array (n_ks_states, n_kpt_this_proc)
    real(dp), contiguous, intent(in) :: ks_energies(:, :)
    !> State occupations array (n_ks_states, n_kpt_this_proc)
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Energy gap value 
    real(dp), intent(out) :: energy_gap

    real(dp), allocatable :: cb_min_all_ranks(:), vb_max_all_ranks(:)
    integer(i32) :: vbm_band_ind, cbm_band_ind, vbm_kpt_ind, cbm_kpt_ind, gap_min_kpt_ind
    
    CALL_ASSERT( all( shape( ks_energies ) == shape( occupations ) ), "Incompatible energies  and occupations arrays provided to get_energy_gap." )
    
    call find_vbm_cbm( 1, size( ks_energies, 1 ), size( ks_energies, 2 ), occupations, &
      ks_energies, vbm_band_ind, cbm_band_ind, vbm_kpt_ind, cbm_kpt_ind, gap_min_kpt_ind )
    allocate(cb_min_all_ranks(mpi_env%procs))
    allocate(vb_max_all_ranks(mpi_env%procs))
    cb_min_all_ranks(mpi_env%rank + 1) = ks_energies(cbm_band_ind, cbm_kpt_ind)
    vb_max_all_ranks(mpi_env%rank + 1) = ks_energies(vbm_band_ind, vbm_kpt_ind)

    call xmpi_allgatherv( mpi_env, cb_min_all_ranks, 1 )
    call xmpi_allgatherv( mpi_env, vb_max_all_ranks, 1 )
    energy_gap = max( 0._dp, minval( cb_min_all_ranks ) - maxval( vb_max_all_ranks ) )
  end subroutine

  !> (private) Wrapper for calling warning
  subroutine wrapper( warning_header, message )
    character(len=*), intent(in) :: warning_header
    character(len=*), intent(in) :: message

    call warning( warning_header // message )
  end subroutine
end module
