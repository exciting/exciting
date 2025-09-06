module rttddft_timings
  use MD, only: MD_timing
  use precision, only: dp
  use asserts, only: assert

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

  !> Type to store timings related to the update of the density in RT-TDDFT
  type, public :: Timing_RTTDDFT_density
    !> timing: total time spent to update the density
    real(dp) :: total
    !> timing: time spent to execute `rhovalk`, `genrhoir`, and eventually 
    !> `mpisumrhoandmag`, see [[update_density]]
    real(dp) :: rho
    !> timing: execution of `symrf`, see [[update_density]]
    real(dp) :: symrf
    !> timing: execution of `rfmtctof`, see [[update_density]]
    real(dp) :: rfmtctof
    !> timing: execution of `addrhocr`, see [[update_density]]
    real(dp) :: addrhocr
    !> timing: execution of `charge`, see [[update_density]]
    real(dp) :: charge
    !> timing: execution of `rhonorm`, see [[update_density]]
    real(dp) :: rhonorm
    !> timing: execution of `from_ks_to_lapw`, see [[update_density]]
    real(dp) :: basis
  contains
    procedure :: reset => reset_Timing_RTTDDFT_density
  end type 

  !> Type to store timings related to the update of the density in RT-TDDFT
  type, public :: Timing_RTTDDFT_potential
    !> timing: total time spent to update the KS potential
    real(dp) :: total
    !> timing: execution of `poteff`, see [[update_potential]]
    real(dp) :: poteff
    !> timing: execution of `genveffig`, see [[update_potential]]
    real(dp) :: genveffig
    !> timing: execution of `genmeffig`, see [[update_potential]]
    real(dp) :: genmeffig
  contains
    procedure :: reset => reset_Timing_RTTDDFT_potential
  end type 

  !> Type to store timings related to the update of the Hamiltonian
  type, public :: Timing_RTTDDFT_hamiltonian
    !> timing: update the hamiltonian
    real(dp) :: total
    !> timing: execution of MT integrals, see [[update_hamiltonian_without_pa_term_lapw]]
    real(dp) :: hmlint
  contains
    procedure :: reset => reset_Timing_RTTDDFT_hamiltonian
  end type 

  type, public :: Timing_RTTDDFT_overlap
    !> timing: update of the overlap via execution of `update_overlap_lapw`, see [[update_overlap_lapw]]
    real(dp) :: total
  contains
    procedure :: reset => reset_Timing_RTTDDFT_overlap
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
    !> timing: evaluation of the time-dependent overlap and berry phase coupling term
    real(dp) :: td_berry
    !> object to store timings spent in the update of the KS hamiltonian
    type(Timing_RTTDDFT_hamiltonian) :: ham
    !> timing: predictor-corrector loop
    real(dp) :: pred_corr
    !> timing: computation of the total energy
    real(dp) :: energy
    !> timing: calculation of the number of excited electrons (per unit cell)
    real(dp) :: n_exc
    !> timing: time for obtaining a screenshot
    real(dp) :: screenshot
    !> timing: time for print operations
    real(dp) :: t_print
  contains
    procedure :: reset => reset_Timing_RTTDDFT
  end type

  !> Type to store timings for Ehrenfest MD
  type, public, extends (MD_timing) :: Timing_Ehrenfest
    !> time to recalculate `pmat` in an Ehrenfest MD step
    real(dp) :: pmat
    !> time to recalculate the hamiltonian matrix in an Ehrenfest MD step
    real(dp) :: ham
    !> time to recalculate the overlap matrix in an Ehrenfest MD step
    type(Timing_RTTDDFT_overlap) :: overlap
  contains
    procedure :: reset => reset_Timing_Ehrenfest
  end type

  type, public :: Timing_RTTDDFT_and_MD
    !> type that contains timings for RT-TDDFT (for evolving KS wavefunctions)
    type(Timing_RTTDDFT) :: t_RTTDDFT
    !> type that contains timings in an Ehrenfest MD
    type(Timing_Ehrenfest) :: t_Ehrenfest
    !> timing: time of each iteration (RT-TDDFT plus MD)
    real(dp) :: t_iteration
  contains
    procedure, public :: reset => reset_Timing_RTTDDFT_and_MD
  end type 


contains
  subroutine set( this, general__, detailed__ )
    class(Print_Timings), intent(inout) :: this
    logical, intent(in) :: general__, detailed__

    if ( detailed__ ) call assert( general__, "detailed timing requested with no general one" )
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

  !> Set every timing value to zero
  pure subroutine reset_Timing_RTTDDFT_and_MD( this )
    class(Timing_RTTDDFT_and_MD), intent(inout) :: this

    call this%t_RTTDDFT%reset()
    call this%t_Ehrenfest%reset()
    this%t_iteration = 0._dp

  end subroutine reset_Timing_RTTDDFT_and_MD

  pure subroutine reset_Timing_RTTDDFT( this )
    class(Timing_RTTDDFT), intent(inout) :: this

    call this%dens%reset()
    call this%pot%reset()
    call this%ham%reset()
    this%wavefunction = 0._dp
    this%current_density = 0._dp
    this%vector_potential = 0._dp
    this%td_berry = 0._dp
    this%pred_corr = 0._dp
    this%energy = 0._dp
    this%n_exc = 0._dp
    this%screenshot = 0._dp
    this%t_print = 0._dp

  end subroutine reset_Timing_RTTDDFT

  pure subroutine reset_Timing_Ehrenfest( this )
    class(Timing_Ehrenfest), intent(inout) :: this

    call this%reset_MD_timing()
    call this%overlap%reset()
    this%pmat = 0._dp
    this%ham = 0._dp

  end subroutine reset_Timing_Ehrenfest

  pure subroutine reset_Timing_RTTDDFT_density( this )
    class(Timing_RTTDDFT_density), intent(inout) :: this
  
    this%total = 0._dp
    this%rho = 0._dp
    this%symrf = 0._dp
    this%rfmtctof = 0._dp
    this%addrhocr = 0._dp
    this%charge = 0._dp
    this%rhonorm = 0._dp
    this%basis = 0._dp
  
  end subroutine reset_Timing_RTTDDFT_density

  pure subroutine reset_Timing_RTTDDFT_potential( this )
    class(Timing_RTTDDFT_potential), intent(inout) :: this

    this%total = 0._dp
    this%poteff = 0._dp
    this%genveffig = 0._dp
    this%genmeffig = 0._dp

  end subroutine reset_Timing_RTTDDFT_potential

  pure subroutine reset_Timing_RTTDDFT_hamiltonian( this )
    class(Timing_RTTDDFT_hamiltonian), intent(inout) :: this

    this%total = 0._dp
    this%hmlint = 0._dp
  end subroutine reset_Timing_RTTDDFT_hamiltonian

  pure subroutine reset_Timing_RTTDDFT_overlap( this )
    class(Timing_RTTDDFT_overlap), intent(inout) :: this

    this%total = 0._dp
  end subroutine

  !> Check the clock (current execution time, in seconds) and store the 
  !> difference between the current time and `ti` (passed as `inout` argument).  
  subroutine timesec_RTTDDFT(ti, duration )
    !> In: initial time in sec
    !> Out: current time in sec
    real(dp),intent(inout) :: ti
    !> Duration, measured as the current time minus `ti`
    real(dp),intent(out) :: duration

    real(dp) :: tf

    call timesec(tf)
    duration = tf - ti
    ti = tf
  end subroutine

end module