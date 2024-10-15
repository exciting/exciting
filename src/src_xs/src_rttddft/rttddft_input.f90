module rttddft_input
  use modinput, only: realTimeTDDFT_type
  use modmpi, only: terminate
  use precision, only: dp, i32
  use propagators, only: propagator_input_elements
  use rttddft_timings, only: Print_Timings
  use rttddft_VectorPotential, only: Vector_Potential

  implicit none

  private

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
    !> Type to encapsulate the attributes of screenshots
    type(screenshot_keys)                   :: screenshots
    !> Type to encapsulate the attributes of pmat
    type(pmat_keys)                         :: pmat
    !> Type to encapsulate the attributes of predictorCorrector
    type(predictorCorrector_keys)           :: predictor_corrector
    !> Type to encapsulate the elements related to the WF propagation
    type(propagator_input_elements)         :: propagator_input
    !> Wether the KS wavefunctions must be normalized in each step 
    logical                                 :: normalize_WF
    !> Print output data every `n_print` steps
    integer(i32)                            :: n_print
    !> Upper limit of time \( t \) - up to which the RT-TDDFT takes place
    real(dp)                                :: t_end
    !> Type that encapsulates if general/detailed information about the RT-TDDFT timings must be printed out
    type(Print_Timings)                     :: printTimings
    !> If `.true.`, calculate of the total energy
    logical                                 :: calculate_total_energy
    !> If `.true.`, calculate of the number of excited electrons
    logical                                 :: calculate_n_exc
    !> If `.true.`, subtract the current density of \(t=0\)
    logical                                 :: subtract_J0
  contains
    procedure :: parse_input => rttddft_input_keys_parse_input
  end type

contains

subroutine rttddft_input_keys_parse_input( this, rt_input, tol, a_vec )
  class(rttddft_input_keys), intent(inout) :: this
  !> Elements and attributes of RT-TDDFT defined in the input file
  type(realTimeTDDFT_type), intent(in) :: rt_input
  !> Tolerance for the methods that need diagonalization
  real(dp), intent(in) :: tol
  !> Type to encapsulate the elements and attributes of laser/vector_potential
  type(Vector_Potential), intent(inout) :: a_vec

  this%normalize_WF = rt_input%normalizeWF
  this%n_print = rt_input%printAfterIterations
  this%t_end = rt_input%endTime
  this%calculate_total_energy = rt_input%calculateTotalEnergy
  this%calculate_n_exc = rt_input%calculateNExcitedElectrons
  this%subtract_J0 = rt_input%subtractJ0
  call this%printTimings%set( rt_input%printTimingGeneral, rt_input%printTimingGeneral .and. rt_input%printTimingDetailed )
  call this%propagator_input%initialize( rt_input%propagator, rt_input%timeStep, rt_input%TaylorOrder, tol )
  call a_vec%initialize( rt_input%laser, rt_input%vectorPotentialSolver )
  
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

end module