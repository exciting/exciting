!> Module for vector vector multiplication.
!> The interfaces combine wrappers for the LAPACK routines
!> ddot, zdotc, zdotu
module vector_multiplication
  use precision, only: dp
  use constants, only: zone
  use asserts, only: assert
  use lapack_f95_interfaces, only: ddot, zdotc, zdotu

  implicit none

  private
  public :: dot_multiply

  !> Spacing between the elements of an vector.
  integer, parameter :: storage_spacing = 1


  !> Calculate the scalar or dot product between two vectors:
  !> \[ 
  !>   \mathbf{a}^dagger \cdot \mathbf{b}
  !> \]
  interface dot_multiply
    module procedure :: dot_multiplication_real_dp, &
                        dot_multiplication_complex_dp, &
                        dot_multiplication_real_complex_dp, &
                        dot_multiplication_complex_real_dp
  end interface dot_multiply


contains


  !> Calculate the scalar or dot product between two real vectors.
  !> \[
  !>   \mathbf{a}^T \cdot \mathbf{b}
  !> \]
  real(dp) function dot_multiplication_real_dp(a, b) 
    !> Input vectors
    real(dp), intent(in) :: a(:), b(:)

    call assert(size(a) == size(b), &
              'a and b must have the same size.')

    dot_multiplication_real_dp = ddot(size(a), a, storage_spacing, b, storage_spacing)
  end function dot_multiplication_real_dp


  !> Calculate the scalar or dot product between two real vectors.
  !> \[
  !>   \mathbf{a}^\dagger \cdot \mathbf{b}
  !> \]
  complex(dp) function dot_multiplication_complex_dp(a, b, conjg_a) 
    !> Input vectors
    complex(dp), intent(in) :: a(:), b(:)
    !> Decide if a is takens as conjungated (default) or tranposed.
    logical, intent(in), optional :: conjg_a

    logical :: conjg_a_ 

    conjg_a_ = .true.
    if (present(conjg_a)) conjg_a_ = conjg_a

    call assert(size(a) == size(b), &
              'a and b must have the same size.')

    if (conjg_a_) then
      dot_multiplication_complex_dp = zdotc(size(a), a, storage_spacing, b, storage_spacing)
    else 
      dot_multiplication_complex_dp = zdotu(size(a), a, storage_spacing, b, storage_spacing)
    end if
  end function dot_multiplication_complex_dp


  !> Calculate the scalar or dot product between two real vectors.
  !> \[
  !>   \mathbf{a}^T \cdot \mathbf{b}
  !> \]
  complex(dp) function dot_multiplication_real_complex_dp(a, b)
    !> Real input vector
    real(dp), intent(in) :: a(:)
    !> Complex input vector
    complex(dp), intent(in) :: b(:)

    call assert(size(a) == size(b), &
              'a and b must have the same size.')

    dot_multiplication_real_complex_dp = zdotu(size(a), zone * a, storage_spacing, b, storage_spacing)
  end function dot_multiplication_real_complex_dp


  !> Calculate the scalar or dot product between two real vectors.
  !> \[
  !>   \mathbf{a}^\dagger \cdot \mathbf{b}
  !> \]
  complex(dp) function dot_multiplication_complex_real_dp(a, b, conjg_a)
    !> Complex input vector
    complex(dp), intent(in) :: a(:)
    !> Real input vector
    real(dp), intent(in) :: b(:)
    !> Decide if a is takens as conjungated (default) or tranposed.
    logical, intent(in), optional :: conjg_a

    logical :: conjg_a_ 

    conjg_a_ = .true.
    if (present(conjg_a)) conjg_a_ = conjg_a

    call assert(size(a) == size(b), &
              'a and b must have the same size.')

    if (conjg_a_) then
      dot_multiplication_complex_real_dp = zdotc(size(a), a, storage_spacing, zone * b, storage_spacing)
    else 
      dot_multiplication_complex_real_dp = zdotu(size(a), a, storage_spacing, zone * b, storage_spacing)
    end if
  end function dot_multiplication_complex_real_dp

end module vector_multiplication