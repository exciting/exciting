module rttddft_timings
  use MD, only: MD_timing
  use precision, only: dp

  implicit none

  private

  public :: timesec_RTTDDFT

  !> Structure to store information about which timings are expected to be printed out
  type, public :: Print_Timings
    private
    !> When `true`, general timings should be printed out
    logical :: general_ = .false.
    !> When `true`, apart from the basic timings, detailed timings should be printed out
    logical :: detailed_ = .false.
  contains
    procedure, public :: set
    procedure, public :: get
    procedure, public :: general
    procedure, public :: detailed
  end type

  !> Type to store timings for Ehrenfest MD
  type, public, extends (MD_timing) :: Timing_Ehrenfest
    !> if `.True.`, it means that an MD step was conducted
    !> This is needed since the time step for MD is a multiple of the time step for RT-TDDFT
    logical  :: MD_was_carried_out
    !> time to recalculate `pmat` in an Ehrenfest MD step
    real(dp) :: pmat
    !> time to recalculate the hamiltonian and overlap matrices in an Ehrenfest MD step
    real(dp) :: hamoverl
  end type

  !> Type to store timings related to the update of the density in RT-TDDFT
  type, public :: Timing_RTTDDFT_density
    !> timing: total time spent to update the density
    real(dp) :: total
    !> timing: time spent to execute `rhovalk`, `genrhoir`, and eventually 
    !> `mpisumrhoandmag`, see [[UpdateDensity]]
    real(dp) :: rho
    !> timing: execution of `symrf`, see [[UpdateDensity]]
    real(dp) :: symrf
    !> timing: execution of `rfmtctof`, see [[UpdateDensity]]
    real(dp) :: rfmtctof
    !> timing: execution of `addrhocr`, see [[UpdateDensity]]
    real(dp) :: addrhocr
    !> timing: execution of `charge`, see [[UpdateDensity]]
    real(dp) :: charge
    !> timing: execution of `rhonorm`, see [[UpdateDensity]]
    real(dp) :: rhonorm
  end type 

  !> Type to store timings related to the update of the density in RT-TDDFT
  type, public :: Timing_RTTDDFT_potential
    !> timing: total time spent to update the KS potential
    real(dp) :: total
    !> timing: execution of `poteff`, see [[uppot]]
    real(dp) :: poteff
    !> timing: execution of `genveffig`, see [[uppot]]
    real(dp) :: genveffig
    !> timing: execution of `genmeffig`, see [[uppot]]
    real(dp) :: genmeffig
  end type 

  !> Type to store timings related to the update of the Hamiltonian
  type, public :: Timing_RTTDDFT_hamiltonian
    !> timing: update the hamiltonian
    real(dp) :: total
    !> timing: execution of `hmlint`, see [[UpdateHam]]
    real(dp) :: hmlint
    !> timing: time spent after executing `hmlint` until the update of the 
    !> hamiltonian has been concluded, see [[UpdateHam]]
    real(dp) :: rest
  end type 

  !> This type stores the time (in seconds) spent in the procedures of RT-TDDFT
  type, public :: Timing_RTTDDFT
    !> timing: evolution of the wavefunction
    real(dp) :: wavefunction
    !> object to store timings spent in the update of the density
    type(Timing_RTTDDFT_density) :: dens
    !> object to store timings spent in the update of the KS potential
    type(Timing_RTTDDFT_potential) :: pot
    !> timing: update of the paramagnetic component of the induced current density
    real(dp) :: current_density
    !> timing: update of the vector potential
    real(dp) :: vector_potential
    !> object to store timings spent in the update of the KS potential
    type(Timing_RTTDDFT_hamiltonian) :: ham
    !> timing: predictor-corrector loop
    real(dp) :: pred_corr
    !> timing: computation of the total energy
    real(dp) :: energy
    !> timing: calculation of the number of excited electrons (per unit cell)
    real(dp) :: n_exc
    !> timing: time for obtaining a screenshot
    real(dp) :: screenshot
  end type

  type, public :: Timing_RTTDDFT_and_MD
    !> type that contains timings for RT-TDDFT (for evolving KS wavefunctions)
    type(Timing_RTTDDFT)   :: t_RTTDDFT
    !> type that contains timings in an Ehrenfest MD
    type(Timing_Ehrenfest) :: t_Ehrenfest
    !> timing: time of each iteration (RT-TDDFT plus MD)
    real(dp) :: t_iteration
  end type 


contains
  pure subroutine set( this, general__, detailed__ )
    class(Print_Timings), intent(inout) :: this
    logical, intent(in) :: general__, detailed__

    this%general_ = general__
    this%detailed_ = detailed__
  end subroutine

  pure subroutine get( this, general__, detailed__ )
    class(Print_Timings), intent(in) :: this
    logical, intent(out) :: general__, detailed__
    
    general__ = this%general_
    detailed__ = this%detailed_
  end subroutine

  pure logical function general( this )
    class(Print_Timings), intent(in) :: this
    general = this%general_
  end function

  pure logical function detailed( this )
    class(Print_Timings), intent(in) :: this
    detailed = this%detailed_
  end function

  !> Check the clock (current execution time, in seconds) and store the 
  !> difference between the current time and `ti` (passed as `inout` argument).  
  subroutine timesec_RTTDDFT(ti, duration )
    !> In: initial time in sec
    !> Out: current time in sec
    real(dp),intent(inout) :: ti
    !> Duration, measured as the current time minus `ti`
    real(dp),intent(out)   :: duration

    real(dp) :: tf

    call timesec(tf)
    duration = tf - ti
    ti = tf
  end subroutine

end module