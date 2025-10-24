!> Module to encapsulate arrays used in RT-TDDFT
! TODO(Ronaldo): Extend to distributed arrays for Scalapack
module rttddft_arrays
  use asserts, only: assert
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
    procedure, public  :: assert_allocated => matrix_set_assert_allocated
    procedure, public  :: allocate_array => matrix_set_allocate_array
    procedure, public  :: copy_object => matrix_set_copy
    procedure, private :: copy_diagonal_real_dp => matrix_set_copy_diagonal_real_dp
    generic, public    :: copy_from => copy_object, copy_diagonal_real_dp
    procedure, public  :: deallocate_if_allocated => matrix_set_deallocate_if_allocated
    procedure, private :: subtract_diagonal_real_dp => matrix_set_subtract_diagonal_real_dp
    procedure, private :: subtract_object => matrix_set_subtract
    generic, public    :: subtract => subtract_object, subtract_diagonal_real_dp
  end type

  !> Type intended to encapsulate generic matrices
  type, public, extends(matrix_set) :: generic_matrix_set
  end type

  !> Type intended to encapsulate hermitian matrices
  type, public, extends(generic_matrix_set) :: hermitian_matrix_set
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

  subroutine matrix_set_assert_allocated( this )
    class(matrix_set), intent(in) :: this

    call assert( allocated(this%array), "array must be allocated")
  end subroutine

  subroutine matrix_set_copy( this, other )
    class(matrix_set), intent(inout) :: this
    class(matrix_set), intent(in) :: other

    call this%assert_allocated()
    call other%assert_allocated()
    this%array = other%array
  end subroutine

  subroutine matrix_set_copy_diagonal_real_dp( this, diagonal )
    class(matrix_set), intent(inout) :: this
    real(dp), intent(in) :: diagonal(:, :)

    integer(i32) :: j, k, ji, jf, ki, kf

    call this%assert_allocated()
    call assert( all([size(this%array, 1), size(this%array, 3) ] == shape(diagonal)), "shape mismatch" )
    call assert( size(this%array, 1) == size(this%array, 2), "not square" )
    ki = lbound(this%array, 3); kf = ubound(this%array, 3)
    ji = lbound(this%array, 1); jf = ubound(this%array, 1) 
    do k = ki, kf
      do j = ji, jf
        this%array(j, j, k) = cmplx( diagonal(j-ji+1, k-ki+1), kind = dp )
      end do
    end do
  end subroutine

  subroutine matrix_set_subtract_diagonal_real_dp( this, diagonal )
    class(matrix_set), intent(inout) :: this
    real(dp), intent(in) :: diagonal(:, :)

    integer(i32) :: j, k, ji, jf, ki, kf

    call this%assert_allocated()
    call assert( all([size(this%array, 1), size(this%array, 3) ] == shape(diagonal)), "shape mismatch" )
    call assert( size(this%array, 1) == size(this%array, 2), "not square" )
    ki = lbound(this%array, 3); kf = ubound(this%array, 3)
    ji = lbound(this%array, 1); jf = ubound(this%array, 1) 
    do k = ki, kf
      do j = ji, jf
        this%array(j, j, k) = this%array(j, j, k) - diagonal(j-ji+1, k-ki+1)
      end do
    end do
  end subroutine

  subroutine matrix_set_subtract( this, other )
    class(matrix_set), intent(inout) :: this
    class(matrix_set), intent(in) :: other

    call this%assert_allocated()
    call other%assert_allocated()
    call assert( all(shape(this%array) == shape(other%array)), "shape mismatch" )
    this%array = this%array - other%array
  end subroutine
end module