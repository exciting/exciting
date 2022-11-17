!> Initialise variables required by the SIRIUS external library 
module sirius_init
  use precision, only: dp
  use iso_c_binding, only: c_ptr, C_NULL_PTR
  use modinput, only: input_type
  use modmpi, only: terminate_mpi_env, mpiglobal

  implicit none
  private

  !> use sirius to generate G+k vector related arrays
  logical, public, parameter :: use_sirius_gkvec = .true.
  !> sirius error message
  character(len=100), public, parameter :: sirius_error = "Not compiled with SIRIUS"

  !> Alias for input%groundstate%sirius
  type, public :: sirius_options_t
    !> use sirius to find the eigen states
    logical :: use_eigen_states
    !> use sirius to construct charge density
    logical :: use_density
    !> use sirius to generate XC related functions and energies
    logical :: use_xc
    !> use sirius to generate Hartree potential
    logical :: use_v_hartree 
    logical :: use_initial_rho
    !> use sirius to compute structure factors
    logical :: use_sfacg
    !> use sirius to compute characteristic function (AKA step function)
    logical :: use_c_function
  contains
    procedure :: initialise => initialise_sirius
    !procedure :: get_use_eigen_states
  end type sirius_options_t

  !TODO(Alex) This should be moved!!!
  type(sirius_options_t), public :: sirius_options

contains

  !> Initialise sirius input options
  subroutine initialise_sirius(this, input)
    !> Instance of sirius options
    class(sirius_options_t), intent(inout) :: this
    !> exciting input 
    type(input_type), intent(in) ::input

    if (.not. associated(input%groundstate%sirius)) then
      call terminate_mpi_env(mpiglobal, "Calling sirius initialisation but no sirius option &
           provided in exciting input ")
    endif

    this%use_eigen_states = input%groundstate%sirius%eigenstates
    this%use_density      = input%groundstate%sirius%density
    this%use_xc           = input%groundstate%sirius%xc
    this%use_v_hartree    = input%groundstate%sirius%vha
    this%use_initial_rho  = input%groundstate%sirius%densityinit
    this%use_sfacg        = input%groundstate%sirius%sfacg
    this%use_c_function   = input%groundstate%sirius%cfun

  end subroutine initialise_sirius

end module sirius_init
