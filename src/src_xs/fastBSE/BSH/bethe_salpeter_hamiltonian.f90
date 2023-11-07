! Module for setting up the ISDF Bethe-Salpeter Hamiltonian (BSH)
module bethe_salpeter_hamiltonian
  use precision, only: dp
  use constants, only: zone

  use interaction_kernel, only: interaction_kernel_type
  use vexc_isdf_kernel, only: vexc_isdf_kernel_type
  use wscr_isdf_kernel, only: wscr_isdf_kernel_type

  private
  public :: bsh_type, vexc_isdf_kernel_type, wscr_isdf_kernel_type

  type bsh_type
    !> transition_energies part or the BSH
    real(dp), allocatable, private :: transition_energies(:)
    !> Exchange interaction kernel
    class(interaction_kernel_type), private, pointer :: vexc => null()
    !> Screened interaction kernel
    class(interaction_kernel_type), private, pointer :: wscr => null()

    integer :: n_transitions

    contains
    
    procedure :: finalize
    generic :: set => set_diagonal, set_isdf_vexc, set_isdf_wscr
    procedure :: set_diagonal, set_isdf_vexc, set_isdf_wscr
    procedure :: times_vector
  end type bsh_type

  contains

  subroutine finalize(this)
    class(bsh_type), intent(out) :: this

    if(allocated(this%transition_energies)) deallocate(this%transition_energies)
    if(associated(this%vexc)) this%vexc => null()
    if(associated(this%wscr)) this%wscr => null()
  end subroutine

  subroutine set_diagonal(this, transition_energies)
    class(bsh_type), intent(inout) :: this
    real(dp), intent(in) :: transition_energies(:)

    this%transition_energies = transition_energies
    this%n_transitions = size(this%transition_energies)
  end subroutine set_diagonal

  subroutine set_isdf_vexc(this, vexc)
    class(bsh_type), intent(inout) :: this
    type(vexc_isdf_kernel_type), intent(in), target :: vexc

    this%vexc => vexc
  end subroutine set_isdf_vexc
  
  subroutine set_isdf_wscr(this, wscr)
    class(bsh_type), intent(inout) :: this
    type(wscr_isdf_kernel_type), intent(in), target :: wscr

    this%wscr => wscr
  end subroutine set_isdf_wscr

  integer function n_transitions(this)
    !> ISDF BSH type
    class(bsh_type) :: this

    n_transitions = size(this%transition_energies)
  end function n_transitions

  !> Apply the ISDF BSH \( \mathbf{H}_{BSH} \) to a vector \( \mathbf{v}_\text{in} \) such that
  !> \[
  !>   \mathbf{v}_\text{out} = \mathbf{H}_{BSH} \cdot \mathbf{v}_\text{in}
  !>     = \mathbf{D} \cdot \mathbf{v}_\text{in}
  !>     + \mathbf{V} \cdot \mathbf{v}_\text{in}
  !>     + \mathbf{W} \cdot \mathbf{v}_\text{in},
  !> \]
  !> where \( \mathbf{D} \) is the transition_energies part of the BSH, \( \mathbf{V} \) the exchange interaction kernel
  !> and \( \mathbf{W} \) the screened interaction kernel, \( \mathbf{V} \cdot \mathbf{v}_\text{in} \) is
  !> calculated with [[V_times_vector]] and  \( \mathbf{W} \cdot \mathbf{v}_\text{in} \) is calculated with
  !> [[WA_times_vector]]
  subroutine times_vector(this, vector_in, vector_out)
    !> ISDF BSH
    class(bsh_type), intent(inout) :: this
    !> Input vector \( \mathbf{v}_\text{in} \)
    complex(dp), intent(in) :: vector_in(:)
    !> Output vector \( \mathbf{v}_\text{out} \)
    complex(dp), intent(out) :: vector_out(:)

    complex(dp), allocatable :: vector_work(:)

    vector_out = this%transition_energies * vector_in

    if (associated(this%vexc)) then
      vector_work = vector_in
      call this%vexc%times_vector(vector_work)
      vector_out = vector_out + 2 * vector_work
    end if

    if (associated(this%wscr)) then
      vector_work = vector_in
      call this%wscr%times_vector(vector_work)
      vector_out = vector_out - vector_work
    end if
    
  end subroutine times_vector

end module bethe_salpeter_hamiltonian