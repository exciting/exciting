!> Module for vector vector multiplication.
!> The interfaces combine wrappers for the LAPACK routines
!> dnrm2, dznrm2, ddot, zdotc, zdotu, dger,  zgerc, zgeru
module vector_multiplication
  use precision, only: dp
  use constants, only: zone, zzero
  use asserts, only: assert
  use lapack_f95_interfaces, only: dnrm2, dznrm2, ddot, zdotc, zdotu, dger,  zgerc, zgeru

  implicit none

  private
  public :: norm, dot_multiply, outer_product

  !> Spacing between the elements of an vector.
  integer, parameter :: storage_spacing = 1
  real(dp), parameter :: prefactor_real_dp = 1.0_dp
  complex(dp), parameter :: prefactor_complex_dp = zone


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


  !> Calculate the eucledian norm of a vector \( \mathbf{v} = (v_1, \cdots, v_n) \)
  !> \[
  !>   \sqrt{\sum_{i=1}^n v_i^2}
  !> \]
  interface norm 
    module procedure :: norm2_real_dp, norm2_complex_dp
  end interface


  !> Calculate the outer product between two vectors with a subroutine call:
  !> \[
  !>    \mathbf{C}_\text{out} = \mathbf{a} \cdot \text{op}(\mathbf{b}) + \mathbf{C}_\text{in}
  !> \]
  !> op means that \( \mathbf \) is either transposed or transpose-conjugated. 
  !> The result is a matrix \( \mathbf{C}_\text{out} \in \mathbb{K}^{n \times m} \),
  !> given elementwise as
  !> \[
  !>    C_{\text{out}; i, j} = a_i \cdot \text{op}(b_j) + C_{\text{in}; i, j}.
  !> \]
  interface outer_product
    module procedure :: outer_product_real_dp, & 
                        outer_product_complex_dp, &
                        outer_product_real_complex_dp, &
                        outer_product_complex_real_dp
  end interface outer_product

contains

! norm 

  !> Calculate the eucledian norm of a vector \( \mathbf{v} = (v_1, \cdots, v_n) \)
  !> \[
  !>   \sqrt{\sum_{i=1}^n v_i^2}.
  !> \]
  real(dp) function norm2_real_dp(v)
    !> Vector for which the norm is to be calculated
    real(dp), intent(in) :: v(:)

    norm2_real_dp = dnrm2(size(v), v, storage_spacing)
  end function


  !> Calculate the eucledian norm of a vector \( \mathbf{v} = (v_1, \cdots, v_n) \)
  !> \[
  !>   \sqrt{\sum_{i=1}^n v_i^2}.
  !> \]
  real(dp) function norm2_complex_dp(v)
    !> Vector for which the norm is to be calculated
    complex(dp), intent(in) :: v(:)

    norm2_complex_dp = dznrm2(size(v), v, storage_spacing)
  end function

! dot product

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

    dot_multiplication_real_complex_dp = zdotu(size(a), cmplx(a, 0.0_dp, kind=dp), storage_spacing, b, storage_spacing)
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
      dot_multiplication_complex_real_dp = zdotc(size(a), a, storage_spacing, cmplx(b, 0.0_dp, kind=dp), storage_spacing)
    else 
      dot_multiplication_complex_real_dp = zdotu(size(a), a, storage_spacing, cmplx(b, 0.0_dp, kind=dp), storage_spacing)
    end if
  end function dot_multiplication_complex_real_dp

! outer product

  !> Calculate the outer product between two real vectors \(\mathbf{a}\) and \(\mathbf{b}\)
  !> \[
  !>   \mathbf{C} = \mathbf{a} \cdot \mathbf{b}^T + \mathbf{C}.
  !> \]
  !> The result is a size(\(\mathbf{a}\)) \(\times\) size(\(\mathbf{b}\)) real matrix \(\mathbf{C}\).
  subroutine outer_product_real_dp(a, b, C)
    !> Input vectors
    real(dp), intent(in), contiguous :: a(:), b(:)
    !> Output matrix
    real(dp), intent(inout), contiguous :: C(:, :)

    call assert(size(C, dim=1) == size(a), 'The number of rows of C must be the same as the number of elements of a.')
    call assert(size(C, dim=2) == size(b), 'The number of rows of B must be the same as the number of elements of b.')

    call dger(size(a), size(b), prefactor_real_dp, a, storage_spacing, b, storage_spacing, C, size(a))
  end subroutine outer_product_real_dp


  !> Calculate the outer product between two complex vectors \(\mathbf{a}\) and \(\mathbf{b}\)
  !> \[
  !>   \mathbf{C} = \mathbf{a} \cdot \text{op}(\mathbf{b}) + \mathbf{C},
  !> \]
  !> where op means that \(\mathbf{b}\) is either transposed or transpose-conjugated.
  !> The result is a size(\(\mathbf{a}\)) \(\times\) size(\(\mathbf{b}\)) complex matrix \(\mathbf{C}\).
  subroutine outer_product_complex_dp(a, b, C, conjg_b)
    !> Input vectors
    complex(dp), intent(in), contiguous :: a(:), b(:)
    !> Output matrix
    complex(dp), intent(inout), contiguous :: C(:, :)
    !> Take the conjugate or transpose of \(mathbf{b}\) (default is .false.).
    logical, intent(in), optional :: conjg_b

    logical :: conjg_b_

    conjg_b_ = .false.
    if (present(conjg_b)) conjg_b_ = conjg_b

    call assert(size(C, dim=1) == size(a), 'The number of rows of C must be the same as the number of elements of a.')
    call assert(size(C, dim=2) == size(b), 'The number of rows of B must be the same as the number of elements of b.')

    if (conjg_b_) then
      call zgerc(size(a), size(b), prefactor_complex_dp, a, storage_spacing, b, storage_spacing, C, size(a))
    else 
      call zgeru(size(a), size(b), prefactor_complex_dp, a, storage_spacing, b, storage_spacing, C, size(a))
    end if
  end subroutine outer_product_complex_dp


  !> Calculate the outer product between a real vector \(\mathbf{a}\) and a complex vector \(\mathbf{b}\)
  !> \[
  !>   \mathbf{C} = \mathbf{a} \cdot \text{op}(\mathbf{b}) + \mathbf{C},
  !> \]
  !> where op means that \(\mathbf{b}\) is either transposed or transpose-conjugated.
  !> The result is a size(\(\mathbf{a}\)) \(\times\) size(\(\mathbf{b}\)) complex matrix \(\mathbf{C}\).
  subroutine outer_product_real_complex_dp(a, b, C, conjg_b)
    !> Input vectors
    real(dp), intent(in), contiguous :: a(:)
    complex(dp), intent(in), contiguous ::  b(:)
    !> Output matrix
    complex(dp), intent(inout), contiguous :: C(:, :)
    !> Take the conjugate or transpose of \(mathbf{b}\) (default os .false.).
    logical, intent(in), optional :: conjg_b

    logical :: conjg_b_

    conjg_b_ = .false.
    if (present(conjg_b)) conjg_b_ = conjg_b
    
    call assert(size(C, dim=1) == size(a), 'The number of rows of C must be the same as the number of elements of a.')
    call assert(size(C, dim=2) == size(b), 'The number of rows of B must be the same as the number of elements of b.')

    if (conjg_b_) then
      call zgerc(size(a), size(b), prefactor_complex_dp, cmplx(a, 0.0_dp, kind=dp), storage_spacing, b, storage_spacing, C, size(a))
    else 
      call zgeru(size(a), size(b), prefactor_complex_dp, cmplx(a, 0.0_dp, kind=dp), storage_spacing, b, storage_spacing, C, size(a))
    end if
  end subroutine outer_product_real_complex_dp


  !> Calculate the outer product between a real vector \(\mathbf{a}\) and a complex vector \(\mathbf{b}\)
  !> \[
  !>   \mathbf{C} = \mathbf{a} \cdot \mathbf{b}^T + \mathbf{C},
  !> \]
  !> The result is a size(\(\mathbf{a}\)) \(\times\) size(\(\mathbf{b}\)) complex matrix \(\mathbf{C}\).
  subroutine outer_product_complex_real_dp(a, b, C)
    !> Input vectors
    complex(dp), intent(in), contiguous :: a(:)
    real(dp), intent(in), contiguous :: b(:)
    !> Output matrix
    complex(dp), intent(inout), contiguous :: C(:, :)

    call assert(size(C, dim=1) == size(a), 'The number of rows of C must be the same as the number of elements of a.')
    call assert(size(C, dim=2) == size(b), 'The number of rows of B must be the same as the number of elements of b.')

    call zgeru(size(a), size(b), prefactor_complex_dp, a, storage_spacing, cmplx(b, 0.0_dp, kind=dp), storage_spacing, C, size(a))
  end subroutine outer_product_complex_real_dp

end module vector_multiplication
