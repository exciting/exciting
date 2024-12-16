module rttddft_input
  use modinput, only: realTimeTDDFT_type, plot3d_type
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
    !> If `.true.`, calculate and print real-time density
    logical :: print_density
    !> Grid data fot 3D density plots
    type(plot3d_type), pointer :: plot3d => null()
  contains
    final :: destructor_screenshot_keys
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

  type :: eeInteraction_keys
    !> Flag that tells if IPA should be invoked
    logical :: ipa
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
    !> Type to encapsulate eeInteraction-related propagation parameters
    type(eeInteraction_keys) :: eeInteraction
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
  if ( this%screenshots%on ) then
    this%screenshots%n_steps = rt_input%screenshots%niter
    this%screenshots%print_density = associated( rt_input%screenshots%density )
    if ( this%screenshots%print_density ) then

      allocate( this%screenshots%plot3d )
      allocate( this%screenshots%plot3d%box )
      allocate( this%screenshots%plot3d%box%origin )
      allocate( this%screenshots%plot3d%box%pointarray(3) )
      allocate( this%screenshots%plot3d%box%pointarray(1)%point )
      allocate( this%screenshots%plot3d%box%pointarray(2)%point )
      allocate( this%screenshots%plot3d%box%pointarray(3)%point )
      this%screenshots%plot3d = rt_input%screenshots%density%plot3d
  
    end if
  else
    this%screenshots%print_density = .false.
  end if

  this%pmat%read_pmat_from_file = rt_input%pmat%readFromFile
  this%pmat%write_pmat_to_file = rt_input%pmat%writeToFile .and. (.not. this%pmat%read_pmat_from_file)
  this%pmat%force_pmat_hermitian = rt_input%pmat%forceHermitian

  this%predictor_corrector%on = associated( rt_input%predictorCorrector )
  if ( this%predictor_corrector%on ) then
    this%predictor_corrector%tol = rt_input%predictorCorrector%tol
    this%predictor_corrector%max_steps = rt_input%predictorCorrector%maxIterations
  end if

  this%eeInteraction%ipa = ( trim( rt_input%eeInteraction ) == "IPA" )


end subroutine

impure elemental subroutine destructor_screenshot_keys( this )
  type(screenshot_keys), intent(inout) :: this

  integer :: i

  if ( associated( this%plot3d ) ) then
    if ( associated( this%plot3d%box ) ) then
      if ( associated( this%plot3d%box%origin ) ) then
        if ( associated( this%plot3d%box%pointarray ) ) then
          do i = 1, 3
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