module rttddft_input
#include "asserts.fpp"
  use modinput, only: deltadensityplot_type, eigenvalues_type, occupations_type, &
    plot3d_type, projectionCoefficients_type, screenshots_type, input_type
  use precision, only: dp, i32
  use propagators, only: propagator_input_elements
  use rttddft_io, only: restart_format, binary, hdf5, file_handler
  use rttddft_timings, only: Print_Timings
  use rttddft_VectorPotential, only: Vector_Potential
  use xhdf5_error_handling, only: abort_if_not_hdf5

  implicit none

  private

  !> Enum with start mode
  enum, bind(C)
    enumerator :: start_mode
    enumerator :: fromscratch, fromfile
  end enum

  !> Enum with basis set
  enum, bind(C)
    enumerator :: basis_set
    enumerator :: lapwlo, ks
  end enum

  !> Enum with field coupling
  enum, bind(C)
    enumerator :: field_coupling
    enumerator :: velocity_gauge, berry_phase
  end enum

  !> Enum with the level of account for interelectronic interaction 
  enum, bind(C)
    enumerator :: ee_interaction
    enumerator :: ipa, tdh, adft
  end enum

  type :: screenshot_eigenvalues_keys
    !> If `.true.`, evaluate the eigenvalues when taking a screenshot
    logical :: on
    !> Number of eigenvalues to be evaluated 
    integer(i32) :: n_eigenvalues
    !> Tolerance for evaluating the eigenvalues 
    real(dp) :: tol
  contains
    procedure, private :: parse_input => screenshot_eigenvalues_keys_parse
  end type

  type :: screenshot_projectionCoefficients_keys
    !> If `.true.`, obtain the projection coefficients when taking a screenshot
    logical :: on
    !> Encapsulate the attribute `printAbsoluteValue` in `projectionCoefficients`
    logical :: print_absolute_value
    !> Encapsulate the attribute `format` in `projectionCoefficients`
    character(len=:), allocatable :: output_format
  contains
    procedure, private :: parse_input => screenshot_projectionCoefficients_keys_parse
  end type

  type :: screenshot_occupations_keys
    !> If `.true.`, evaluate occupation numbers when taking a screenshot
    logical :: on
    !> Encapsulate the attribute `format` in `occupations`
    character(len=:), allocatable :: output_format
    !> Encapsulate the attribute `textFormat`
    logical :: output_text_format
    !> Encapsulate the attribute `binaryFormat`
    logical :: output_binary_format
  contains
    procedure, private :: parse_input => screenshot_occupations_keys_parse
  end type

  type :: screenshot_density_keys
    !> If `.true.`, print out the electron density when taking a screenshot
    logical :: on
    !> Grid data fot 3D density plots
    type(plot3d_type), pointer :: plot3d => null()
  contains
    procedure, private :: parse_input => screenshot_density_keys_parse
    final :: destructor_screenshot_density_keys
  end type

  type, public :: screenshot_keys
    !> If `.true.`, take screenshots during the RT-TDDFT evolution
    logical :: on
    !> Take a screenshot every `n_steps` number of steps
    integer(i32) :: n_steps
    !> Type to encapsulate the attributes of `eigenvalues` inside an `screenshots` element
    type(screenshot_eigenvalues_keys) :: eigenvalues
    !> Type to encapsulate the attributes of `projectionCoefficients` inside an `screenshots` element
    type(screenshot_projectionCoefficients_keys) :: projection_coefficients
    !> Type to encapsulate the attributes of `occupations` inside an `screenshots` element
    type(screenshot_occupations_keys) :: occupations
    !> Type to encapsulate the attributes of `density` inside an `screenshots` element
    type(screenshot_density_keys) :: density
  contains
    procedure, private :: parse_input => screenshot_input_keys_parse_input
  end type

  type, public :: pmat_keys
    logical :: read_pmat_from_file
    logical :: write_pmat_to_file
    logical :: force_pmat_hermitian
  end type

  type :: predictorCorrector_keys
    !> Maximum number of steps for the predictor-corrector loop
    integer(i32) :: max_steps
    !> Flag that tells if the predictor-corrector scheme is required
    logical :: on
    !> tolerance to escape the predictor-corrector loop
    real(dp) :: tol
  end type

  type :: eeInteraction_keys
    !> Flag that tells if the adiabatic XC part of the effective should be evolved
    logical, private :: evolve_adiabatic_xc
    !> Flag that tells if effective potential w/o XC part should be evolved
    logical, private :: evolve_coulomb
    contains
      procedure :: use_ipa => eeInteraction_use_ipa, &
                   coulomb_only => eeInteraction_coulomb_only
      procedure, private :: eeInteraction_init_from_string
  end type

  !> Type to encapsulate the elements and attributes defined in the input file
  type, public :: rttddft_input_keys
    !> Type to encapsulate the attributes of screenshots
    type(screenshot_keys) :: screenshots
    !> Type to encapsulate the attributes of pmat
    type(pmat_keys) :: pmat
    !> Type to encapsulate the attributes of predictorCorrector
    type(predictorCorrector_keys) :: predictor_corrector
    !> Type to encapsulate the elements related to the WF propagation
    type(propagator_input_elements) :: propagator_input
    !> Type to encapsulate eeInteraction-related propagation parameters
    type(eeInteraction_keys) :: eeInteraction
    !> Whether the KS wavefunctions must be normalized in each step 
    logical :: normalize_WF
    !> Number of low-lying states which will not be evolved
    integer(i32) :: n_frozen
    !> Whether the active KS wavefunctions must be orthogonalized against frozen in each step 
    logical :: orthogonalize_against_frozen
    !> Print output data every `n_print` steps
    integer(i32) :: n_print
    !> Radial step length (used to update the electron density)
    integer(i32) :: l_rad_step
    !> Upper limit of time \( t \) - up to which the RT-TDDFT takes place
    real(dp) :: t_end
    !> Type that encapsulates if general/detailed information about the RT-TDDFT timings must be printed out
    type(Print_Timings) :: printTimings
    !> If `.true.`, calculate of the total energy
    logical :: calculate_total_energy
    !> If `.true.`, calculate of the number of excited electrons
    logical :: calculate_n_exc
    !> If `.true.`, subtract the current density of \(t=0\)
    logical :: subtract_J0
    !> Energy shift  \( \Delta E \) for the scissor operator
    real(dp) :: scissor_shift
    !> If `.true.`, write a restart file every `n_print` steps
    logical, private :: save_state
    !> If `.true.`, update the SOC term
    logical, private :: updateSOC
    !> Identify which basis set will be used for the propagation (see [[basis_set]])
    integer(kind( basis_set )), private :: basis_set
    !> Identify which operator will be used for the coupling with external field (see [[field_coupling]])
    integer(kind( field_coupling )), private :: field_coupling
    !> Identify if which start mode is desired (see [[start_mode]])
    integer(kind( start_mode )), private :: start_mode
    !> Format handler of the checkpoint (restart) files
    type(file_handler) :: restart_file_handler
    !> When restarting a calculation, append this string as extension to the common output files
    !> E.g. if `restart_extension=".SAVE"`, the vector potential is read from `VECTOR_POTENTIAL.OUT.SAVE`
    character(len=:), allocatable :: restart_extension
  contains
    procedure :: parse_input => rttddft_input_keys_parse_input
    procedure :: write_restart => rttddft_input_keys_write_restart
    procedure :: restart_previous_calculation => rttddft_input_keys_restart_previous_calculation
    procedure :: do_from_scratch => rttddft_input_keys_do_from_scratch
    procedure :: update_SOC => rttddft_input_keys_update_SOC
    procedure :: use_ks_basis => rttddft_input_keys_use_ks_basis
    procedure :: use_lapwlo_basis => rttddft_input_keys_use_lapwlo_basis
    procedure :: use_velocity_gauge => rttddft_input_keys_use_velocity_gauge
    procedure :: use_berry_phase => rttddft_input_keys_use_berry_phase
  end type

contains

subroutine rttddft_input_keys_parse_input( this, inp, tol, a_vec )
  class(rttddft_input_keys), intent(inout) :: this
  !> Elements and attributes defined in the input file
  type(input_type), intent(in) :: inp
  !> Tolerance for the methods that need diagonalization
  real(dp), intent(in) :: tol
  !> Type to encapsulate the elements and attributes of laser/vector_potential
  type(Vector_Potential), intent(inout) :: a_vec

  associate( rt_input => inp%xs%realTimeTDDFT )
    this%normalize_WF = rt_input%normalizeWF
    this%orthogonalize_against_frozen = rt_input%orthogonalizeAgainstFrozen
    this%n_frozen = rt_input%numberOfFrozenStates
    this%n_print = rt_input%printAfterIterations
    this%t_end = rt_input%endTime
    this%calculate_total_energy = rt_input%calculateTotalEnergy
    this%calculate_n_exc = rt_input%calculateNExcitedElectrons
    this%subtract_J0 = rt_input%subtractJ0
    call this%printTimings%set( rt_input%printTimingGeneral, rt_input%printTimingGeneral .and. rt_input%printTimingDetailed )
    call this%propagator_input%initialize( rt_input%propagator, rt_input%timeStep, rt_input%TaylorOrder, tol, rt_input%nEigenvectorsEH )
    call a_vec%initialize( rt_input%laser, rt_input%vectorPotentialSolver )
    
    this%screenshots%on = associated( rt_input%screenshots )
    if( this%screenshots%on ) call this%screenshots%parse_input( rt_input%screenshots )

    if ( associated( rt_input%pmat ) ) then
      this%pmat%read_pmat_from_file = rt_input%pmat%readFromFile
      this%pmat%write_pmat_to_file = rt_input%pmat%writeToFile .and. (.not. this%pmat%read_pmat_from_file)
      this%pmat%force_pmat_hermitian = rt_input%pmat%forceHermitian
    end if

    this%predictor_corrector%on = associated( rt_input%predictorCorrector )
    if ( this%predictor_corrector%on ) then
      this%predictor_corrector%tol = rt_input%predictorCorrector%tol
      this%predictor_corrector%max_steps = rt_input%predictorCorrector%maxIterations
    end if

    call this%eeInteraction%eeInteraction_init_from_string( rt_input%eeInteraction )
    this%save_state = rt_input%saveState
    this%basis_set = string_to_basis_set( rt_input%basis )
    this%field_coupling = string_to_field_coupling( rt_input%fieldCoupling )
    this%start_mode = string_to_start_mode( rt_input%do )
    this%restart_file_handler%file_format = string_to_restart_format( rt_input%restartFilesFormat )
    this%restart_extension = trim( rt_input%restartExtension )
  end associate
  if( this%restart_file_handler%file_format == hdf5 ) call abort_if_not_hdf5( &
    message="exciting needs to be compiled with HDF5 to use the RT-TDDFT restart feature with HDF5" )
  this%restart_file_handler%file_name = trim( inp%xs%h5fname )
  this%restart_file_handler%path = trim( inp%xs%h5gname )
  this%l_rad_step = inp%groundstate%lradstep
  ! updateSOC must be false for spin-unpolarized calculations and when SOC is not used
  this%updateSOC = .false.
  if( associated( inp%groundstate%spin ) ) then
    this%updateSOC = inp%groundstate%spin%spinorb .and. inp%xs%realTimeTDDFT%spinPropagation%updateSOC
  end if
  this%scissor_shift = inp%xs%scissor
end subroutine

!> Check whether the velocity gauge will be used for the coupling with external field
pure logical function rttddft_input_keys_use_velocity_gauge( this ) result( check )
  class(rttddft_input_keys), intent(in) :: this
  check = ( this%field_coupling == velocity_gauge )
end function

!> Check whether the dynamical Berry phase approach will be used for the coupling with external field
pure logical function rttddft_input_keys_use_berry_phase( this ) result( check )
  class(rttddft_input_keys), intent(in) :: this
  check = ( this%field_coupling == berry_phase )
end function

!> Check whether SOC should be updated
pure logical function rttddft_input_keys_update_SOC( this ) result( check )
  class(rttddft_input_keys), intent(in) :: this
  check = this%updateSOC
end function

!> (private) Given a string, get the corresponding [[field_coupling]]
function string_to_field_coupling(string) result(r)
  !> String containing the start mode name
  character(len=*), intent(in) :: string
  integer(kind( field_coupling )) :: r

  select case ( trim( string ) )
    case ("velocityGauge")
      r = velocity_gauge
    case ("berryPhase")
      r = berry_phase
    case default
      CALL_ASSERT( .false., "Unrecognized field_coupling")
  end select
end function

!> (private) Given a string, get the corresponding [[ee_interaction]]
function string_to_ee_interaction( string ) result( r )
  !> String containing the approximation name
  character(len=*), intent(in) :: string
  integer(kind( ee_interaction )) :: r

  select case ( trim( string ) )
    case ("IPA")
      r = ipa
    case ("aDFT")
      r = adft
    case ("tdH")
      r = tdh
    case default
      CALL_ASSERT( .false., "Unrecognized ee interaction")
  end select
end function

!> Check whether IPA is used
pure logical function eeInteraction_use_ipa( this ) result( check )
  class(eeInteraction_keys), intent(in) :: this
  check = .not. ( this%evolve_coulomb .or. this%evolve_adiabatic_xc )
end function

!> Check whether only Coulomb potential should be evolved
pure logical function eeInteraction_coulomb_only( this ) result( check )
  class(eeInteraction_keys), intent(in) :: this
  check = this%evolve_coulomb .and. (.not. this%evolve_adiabatic_xc) 
end function

!> Initialize [[eeInteraction_keys]] flags from string 
subroutine eeInteraction_init_from_string( this, string )
  class(eeInteraction_keys), intent(inout) :: this
  !> String containing the approximation name
  character(len=*), intent(in) :: string

  select case( string_to_ee_interaction( string ) )
  case(ipa)
    this%evolve_adiabatic_xc = .false.
    this%evolve_coulomb = .false.
  case(adft)
    this%evolve_adiabatic_xc = .true.
    this%evolve_coulomb = .true.
  case(tdh)
    this%evolve_adiabatic_xc = .false.
    this%evolve_coulomb = .true.
  case default
    CALL_ASSERT( .false., "Unrecognized ee interaction")
  end select
end subroutine

!> Check whether the ks basis will be used for time propagation
pure logical function rttddft_input_keys_use_ks_basis(this) result(check)
  class(rttddft_input_keys), intent(in) :: this
  check = ( this%basis_set == ks )
end function

!> Check whether the LAPW+lo basis will be used for time propagation
pure logical function rttddft_input_keys_use_lapwlo_basis(this) result(check)
  class(rttddft_input_keys), intent(in) :: this
  check = ( this%basis_set == lapwlo )
end function

!> (private) Given a string, get the corresponding [[start_mode]]
function string_to_basis_set(string) result(r)
  !> String containing the start mode name
  character(len=*), intent(in) :: string
  integer(kind( basis_set )) :: r

  select case ( trim( string ) )
    case ("LAPWlo")
      r = lapwlo
    case ("unperturbedKS")
      r = ks
    case default
      CALL_ASSERT( .false., "Unrecognized basis set")
  end select
end function

!> Check whether the RT-TDDFT calculation will run from file
pure logical function rttddft_input_keys_restart_previous_calculation(this) result(check)
  class(rttddft_input_keys), intent(in) :: this
  check = ( this%start_mode == fromfile )
end function

!> Check whether the RT-TDDFT calculation will run from scratch
pure logical function rttddft_input_keys_do_from_scratch(this) result(check)
  class(rttddft_input_keys), intent(in) :: this
  check = ( this%start_mode == fromscratch )
end function

!> (private) Given a string, get the corresponding [[start_mode]]
function string_to_start_mode(string) result(r)
  !> String containing the start mode name
  character(len=*), intent(in) :: string
  integer(kind(start_mode)) :: r

  select case ( trim(string) )
    case ("fromscratch")
      r = fromscratch
    case ("fromfile")
      r = fromfile
    case default
      CALL_ASSERT( .false., "Unrecognized string")
  end select
end function

!> (private) Given a string, get the corresponding [[restart_format]]
function string_to_restart_format(string) result(r)
  !> String containing the start mode name
  character(len=*), intent(in) :: string
  integer(kind( restart_format )) :: r

  select case ( trim(string) )
    case ("binary")
      r = binary
    case ("hdf5")
      r = hdf5
    case default
      CALL_ASSERT( .false., "Unrecognized string")
  end select
end function

!> Returns `.true.`, if a restart output must be written during the RT-TDDFT evolution
pure logical function rttddft_input_keys_write_restart(this) result(r)
  class(rttddft_input_keys), intent(in) :: this
  r = this%save_state
end function

!> Parse the input keys defined in the `screenshots` element
subroutine screenshot_input_keys_parse_input( this, screenshots_input )
  class(screenshot_keys), intent(inout) :: this
  !> Elements and attributes of `screenshots` defined in the input file
  type(screenshots_type), intent(in) :: screenshots_input

  this%n_steps = screenshots_input%niter
  
  this%eigenvalues%on = associated( screenshots_input%eigenvalues )
  if( this%eigenvalues%on ) call this%eigenvalues%parse_input( screenshots_input%eigenvalues )
  
  this%projection_coefficients%on = associated( screenshots_input%projectionCoefficients )
  if( this%projection_coefficients%on ) call this%projection_coefficients%parse_input( screenshots_input%projectionCoefficients )
  
  this%occupations%on = associated( screenshots_input%occupations )
  if( this%occupations%on ) call this%occupations%parse_input( screenshots_input%occupations )

  this%density%on = associated( screenshots_input%deltadensityplot )
  if( this%density%on ) call this%density%parse_input( screenshots_input%deltadensityplot )

  ! Turn off screenshots if no property is required
  if( .not. ( this%eigenvalues%on .or. this%projection_coefficients%on .or. this%occupations%on .or. this%density%on ) ) this%on = .false.
end subroutine


!> Parse the input keys defined in the `eigenvalues` element
pure subroutine screenshot_eigenvalues_keys_parse( this, eigenvalues_input )
  class(screenshot_eigenvalues_keys), intent(inout) :: this
  !> Elements and attributes of `eigenvalues` defined in the input file
  type(eigenvalues_type), intent(in) :: eigenvalues_input

  this%n_eigenvalues = eigenvalues_input%nEigenvalues
  this%tol = eigenvalues_input%tolerance
end subroutine


!> Parse the input keys defined in the `projectionCoefficients` element
pure subroutine screenshot_projectionCoefficients_keys_parse( this, projectionCoefficients_input )
  class(screenshot_projectionCoefficients_keys), intent(inout) :: this
  !> Elements and attributes of `eigenvalues` defined in the input file
  type(projectionCoefficients_type), intent(in) :: projectionCoefficients_input

  this%print_absolute_value = projectionCoefficients_input%printAbsoluteValue
  this%output_format = projectionCoefficients_input%format
end subroutine


!> Parse the input keys defined in the `occupations` element
pure subroutine screenshot_occupations_keys_parse( this, occupations_input )
  class(screenshot_occupations_keys), intent(inout) :: this
  !> Elements and attributes of `occupations` defined in the input file
  type(occupations_type), intent(in) :: occupations_input

  this%output_format = occupations_input%format
  this%output_text_format = occupations_input%textFormat
  this%output_binary_format = occupations_input%binaryFormat
  ! turn off occupations if both attributes are false
  if( (.not. this%output_text_format) .and. (.not. this%output_binary_format) ) this%on = .false.
end subroutine


!> Parse the input keys defined in the `deltadensityplot` element
subroutine screenshot_density_keys_parse( this, density_input )
  class(screenshot_density_keys), intent(inout) :: this
  !> Elements and attributes of `deltadensityplot` defined in the input file
  type(deltadensityplot_type), intent(in) :: density_input

  allocate( this%plot3d )
  allocate( this%plot3d%box )
  allocate( this%plot3d%box%origin )
  allocate( this%plot3d%box%pointarray(3) )
  allocate( this%plot3d%box%pointarray(1)%point )
  allocate( this%plot3d%box%pointarray(2)%point )
  allocate( this%plot3d%box%pointarray(3)%point )
  this%plot3d = density_input%plot3d
end subroutine


impure elemental subroutine destructor_screenshot_density_keys( this )
  type(screenshot_density_keys), intent(inout) :: this

  integer(i32) :: i

  if ( associated( this%plot3d ) ) then
    if ( associated( this%plot3d%box ) ) then
      if ( associated( this%plot3d%box%origin ) ) then
        if ( associated( this%plot3d%box%pointarray ) ) then
          do i = 1, size(this%plot3d%box%pointarray)
            if ( associated( this%plot3d%box%pointarray(i)%point ) ) &
              deallocate( this%plot3d%box%pointarray(i)%point )
          end do
          deallocate( this%plot3d%box%pointarray )
        end if
        deallocate( this%plot3d%box%origin )
      end if
      deallocate( this%plot3d%box )
    end if
    deallocate( this%plot3d )
  end if
end subroutine 

end module