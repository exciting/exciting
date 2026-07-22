!> Module with the generic vector field type
module rttddft_VectorField
#include "asserts.fpp"
  use precision, only: dp

  implicit none

  private

  public :: cartesian_direction, direction, x, y, z

  !> Generic type for vector fields
  type, abstract, public :: Uniform_Vector_Field
    !> `x`, `y` and `z` components
    real(dp) :: components(3) = 0._dp
  contains
    procedure, public :: add_vector
  end type

  !> Enum with the cartesian directions
  enum, bind(C)
    enumerator :: direction
    enumerator :: x=1, y=2, z=3
  end enum

contains
  !> Return an enum with the cartesian direction
  function cartesian_direction( char ) result( this )
    !> char to be converted to enum
    character, intent(in) :: char
    !> resulting enum
    integer(kind(direction)) :: this

    CALL_ASSERT( char=='x' .or. char=='y' .or. char=='z', 'Invalid direction')

    select case( char )
      case('x')
        this = x
      case('y')
        this = y
      case('z')
        this = z
    end select
  end function

  !> Add a vector to `this%components`
  pure subroutine add_vector( this, vector )
    class(Uniform_Vector_Field), intent(inout) :: this
    !> Vector to be added
    real(dp), intent(in) :: vector(3)

    this%components = this%components + vector
  end subroutine

end module