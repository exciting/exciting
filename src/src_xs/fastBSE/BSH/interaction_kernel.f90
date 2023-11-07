! Created by  on 11/03/2022.

module interaction_kernel
  use precision, only: dp
  use bse_transitions, only: transition_type

  type, abstract :: interaction_kernel_type
    !> Sizes
    integer :: n_transitions, n_k
    type(transition_type) :: transitions
    complex(dp), allocatable :: fill_zeros(:)

  contains
    procedure(finalize_kernel), deferred :: finalize
    procedure(kernel_times_vector), deferred  :: times_vector
  end type interaction_kernel_type

  abstract interface
    subroutine finalize_kernel(this)
      import interaction_kernel_type

      class(interaction_kernel_type), intent(inout) :: this
    end subroutine finalize_kernel

    subroutine kernel_times_vector(this, vector)
      import dp
      import interaction_kernel_type

      class(interaction_kernel_type), target, intent(inout) :: this
      complex(dp), contiguous, target, intent(inout) :: vector(:)
    end subroutine kernel_times_vector
  end interface

  contains

  

end module interaction_kernel