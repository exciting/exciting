!> Module to encapsulate arrays used in RT-TDDFT
! TODO(Ronaldo): Extend to distributed arrays for Scalapack
module rttddft_arrays
  use constants, only: zzero
  use precision, only: dp, i32

  implicit none

  private

  !> Abstract type that should be extended by any concrete matrix set
  type, public, abstract :: matrix_set
    private
    complex(dp), public, allocatable :: array(:, :, :)
  contains
    private
    procedure, public :: allocate_array => matrix_set_allocate_array
    procedure, public :: deallocate_if_allocated => matrix_set_deallocate_if_allocated
  end type

  !> Type intended to encapsulate hermitian matrices
  type, public, extends(matrix_set) :: hermitian_matrix_set
  end type

  !> Type intended to encapsulate positive definite matrices
  type, public, extends(hermitian_matrix_set) :: positive_matrix_set
  end type

  !> Type intended to encapsulate identity matrices
  type, public, extends(positive_matrix_set) :: identity_matrix_set
  end type

contains
  pure subroutine matrix_set_allocate_array( this, lbounds, ubounds )
    class(matrix_set), intent(inout) :: this
    integer(i32), intent(in) :: lbounds(3)
    integer(i32), intent(in) :: ubounds(3)

    call this%deallocate_if_allocated()
    allocate( this%array(lbounds(1):ubounds(1), &
      lbounds(2):ubounds(2), lbounds(3):ubounds(3) ), source=zzero )
  end subroutine

  pure subroutine matrix_set_deallocate_if_allocated( this )
    class(matrix_set), intent(inout) :: this
    
    if( allocated(this%array) ) deallocate( this%array )
  end subroutine
end module