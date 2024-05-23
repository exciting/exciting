module rttddft_input
  use modinput, only: realTimeTDDFT_type
  use modmpi, only: terminate
  use precision, only: dp, i32
  use rttddft_Wavefunction, only: propagator_types, propagator_type, propagator_keys
  use rttddft_VectorPotential, only: solver_types, solver_type, euler

  implicit none

  private

  character(len=*), parameter :: field_total = 'total'
  character(len=*), parameter :: field_external = 'external'

  !> Enum with the type of applied field
  !> There are 2 possibilities that the applied field can assume: "total" or "external"
  enum, bind(C)
    enumerator :: applied_field
    enumerator :: total, external
  end enum

  !> Type to encapsulate the elements and attributes defined in the laser element
  type :: laser_keys
    !> Type of field used for the vector potential
    integer(kind(applied_field)), private   :: field_type
  contains
    procedure         :: is_field_type_external => laser_is_field_type_external
    procedure         :: is_field_type_total => laser_is_field_type_total
  end type

  type :: screenshot_keys
    !> If `.true.`, take screenshots during the RT-TDDFT evolution
    logical :: on
    !> Take a screenshot every `n_steps` number of steps
    integer(i32) :: n_steps
  end type

  type, public :: pmat_keys
    logical :: read_pmat_from_file
    logical :: write_pmat_to_file
    logical :: force_pmat_hermitian
  end type

  type :: predictorCorrector_keys
    !> Maximum number of steps for the predictor-corrector loop
    integer                   :: max_steps
    !> Flag that tells if the predictor-corrector scheme is required
    logical                   :: on
    !> tolerance to escape the predictor-corrector loop
    real(dp)                  :: tol
  end type

  !> Type to encapsulate the elements and attributes defined in the input file
  type, public :: rttddft_input_keys
    !> Type to encapsulate the elements and attributes of laser
    type(laser_keys)                        :: laser
    !> Type to encapsulate the attributes of screenshots
    type(screenshot_keys)                   :: screenshots
    !> Type to encapsulate the attributes of pmat
    type(pmat_keys)                         :: pmat
    !> Type to encapsulate the attributes of predictorCorrector
    type(predictorCorrector_keys)           :: predictor_corrector
    !> Type to encapsulate the elements related to the WF propagation
    type(propagator_keys)                   :: propagator
    !> Print output data every `n_print` steps
    integer(i32)                            :: n_print
    !> Upper limit of time \( t \) - up to which the RT-TDDFT takes place
    real(dp)                                :: t_end
    !> Type of solver used for the vector potential
    integer(kind(solver_types))             :: vector_potential_solver
    !> If `.true.`, print out general information about the RT-TDDFT timings
    logical                                 :: timings_general
    !> If `.true.`, print out detailed information about the RT-TDDFT timings
    logical                                 :: timings_detailed
    !> If `.true.`, calculate of the total energy
    logical                                 :: calculate_total_energy
    !> If `.true.`, calculate of the number of excited electrons
    logical                                 :: calculate_n_exc
    !> If `.true.`, subtract the current density of \(t=0\)
    logical                                 :: subtract_J0
  contains
    procedure         :: parse_input => rttddft_input_keys_parse_input
    procedure         :: is_field_type_external => rttddft_input_is_field_type_external
    procedure         :: is_field_type_total => rttddft_input_is_field_type_total
    procedure         :: is_solver_euler => rttddft_input_is_solver_euler
  end type

contains

subroutine rttddft_input_keys_parse_input( this, rt_input, tol )
  class(rttddft_input_keys), intent(inout) :: this
  !> Elements and attributes of RT-TDDFT defined in the input file
  type(realTimeTDDFT_type), intent(in) :: rt_input
  !> Tolerance for the methods that need diagonalization
  real(dp), intent(in) :: tol

  this%n_print = rt_input%printAfterIterations
  this%t_end = rt_input%endTime
  this%vector_potential_solver = solver_type( rt_input%vectorPotentialSolver )
  this%calculate_total_energy = rt_input%calculateTotalEnergy
  this%calculate_n_exc = rt_input%calculateNExcitedElectrons
  this%subtract_J0 = rt_input%subtractJ0
  this%timings_general = rt_input%printTimingGeneral
  this%timings_detailed = this%timings_general .and. rt_input%printTimingDetailed

  this%propagator%name = propagator_type( rt_input%propagator )
  this%propagator%time_step = rt_input%timeStep
  this%propagator%normalize_WF = rt_input%normalizeWF
  this%propagator%order_taylor = rt_input%TaylorOrder
  this%propagator%tol = tol

  this%laser%field_type = field_type( rt_input%laser%fieldType )
  
  this%screenshots%on = associated( rt_input%screenshots )
  if( this%screenshots%on ) this%screenshots%n_steps = rt_input%screenshots%niter

  this%pmat%read_pmat_from_file = rt_input%pmat%readFromFile
  this%pmat%write_pmat_to_file = rt_input%pmat%writeToFile .and. (.not. this%pmat%read_pmat_from_file)
  this%pmat%force_pmat_hermitian = rt_input%pmat%forceHermitian

  this%predictor_corrector%on = associated( rt_input%predictorCorrector )
  if ( this%predictor_corrector%on ) then
    this%predictor_corrector%tol = rt_input%predictorCorrector%tol
    this%predictor_corrector%max_steps = rt_input%predictorCorrector%maxIterations
  end if

end subroutine


pure logical function laser_is_field_type_external( this )
  class(laser_keys), intent(in) :: this

  laser_is_field_type_external = ( this%field_type == external )
end function


pure logical function laser_is_field_type_total( this )
  class(laser_keys), intent(in) :: this

  laser_is_field_type_total = ( this%field_type == total )
end function


pure logical function rttddft_input_is_field_type_external( this )
  class(rttddft_input_keys), intent(in) :: this

  rttddft_input_is_field_type_external = this%laser%is_field_type_external()
end function


pure logical function rttddft_input_is_field_type_total( this )
  class(rttddft_input_keys), intent(in) :: this

  rttddft_input_is_field_type_total = this%laser%is_field_type_total()
end function


pure logical function rttddft_input_is_solver_euler( this )
  class(rttddft_input_keys), intent(in) :: this

  rttddft_input_is_solver_euler = ( this%vector_potential_solver == euler )
end function

function field_type( name ) result( field )
  character(len=*), intent(in) :: name
  integer(kind(applied_field)) :: field

  select case( trim( name ) )
    case( field_total )
      field = total
    case( field_external )
      field = external
    case default
      call terminate('unknow laser field type')
  end select
end function


end module