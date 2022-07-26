module inverse
  use precision, only: dp
  use asserts, only: assert
  use math_utils, only: is_square
  use lu_factorization, only: xgetrf2, xgetri
  use determinant, only: determinant_LU

  private
  public :: invert_LU, inverse_LU

  !> Default value of the tolerance for the determinant taken as zero.
  real(dp), parameter :: default_det_tol = 1.e-8_dp

  !> Invert the matrix \( \mathbf{A} \) by calculating the LU factorization with **[[lu_row_pivot]]**,
  !> inverting \( \mathbf{U} \) and solving the system of linear equations
  !> \[
  !>    \mathbf{A}^{-1} \cdot \mathbf{L} = \mathbf{U}^{-1}
  !> \]
  !> Use this interface if \(\mathbf{A}\) can be overwritten.
  interface invert_LU
    module procedure invert_LU_real_dp, &
                     invert_LU_complex_dp
  end interface invert_LU

  !> Invert the matrix \( \mathbf{A} \) by calculating the LU factorization with **[[lu_row_pivot]]**,
  !> inverting \( \mathbf{U} \) and solving the system of linear equations
  !> \[
  !>    \mathbf{A}^{-1} \cdot \mathbf{L} = \mathbf{U}^{-1}
  !> \]
  !> Use this interface if \(\mathbf{A}\) should be unchanged.
  interface inverse_LU
    module procedure inverse_LU_real_dp, &
                     inverse_LU_complex_dp
  end interface inverse_LU

  contains

  !---------------------------------------------------------------------------------------------------------------------
  ! invert_LU, inverse_LU

  !> Invert the matrix \( \mathbf{A} \) with **[[dgetrf2_wrapper]]** and **[[dgetri_wrapper]]**.
  !> If this interface is used, the array `A` contains afterwards the matrix \( \mathbf{A}^{-1} \).
  !>
  !> This is a wrapper for the LAPACK routines dgetrf2 and dgetri.
  subroutine invert_LU_real_dp(A)
    !> On entry, matrix to invert \( \mathbf{A} \) and
    !> on exit inverted matrix \(mathbf{A}^{-1}\)
    real(dp), intent(inout), contiguous :: A(:, :)

    integer, allocatable :: row_pivot_map(:)

    call assert(is_square(A), 'A is not a square matrix.')
    call assert(abs(determinant_LU(A)) > default_det_tol, 'Determinant of A is zero.')

    allocate(row_pivot_map(size(A, dim=1)))

    call xgetrf2(A, row_pivot_map)
    call xgetri(A, row_pivot_map)
  end subroutine invert_LU_real_dp


  !> See **[[invert_LU_real_dp]]**. If this interface is used, the array `A` is unchanged and the matrix
  !> \( \mathbf{A}^{-1} \) is copied to a new array.
  !>
  !> This is a wrapper for the LAPACK routines dgetrf2 and dgetri.
  function inverse_LU_real_dp(A) result(A_inv)
    !> Matrix to invert \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)

    real(dp), allocatable :: A_inv(:, :)

    A_inv = A
    call invert_LU_real_dp(A_inv)
  end function inverse_LU_real_dp


  !> Invert the matrix \( \mathbf{A} \) with **[[zgetrf2_wrapper]]** and **[[zgetri_wrapper]]**.
  !> If this interface is used, the array `A` contains afterwards the matrix \( \mathbf{A}^{-1} \).
  !>
  !> This is a wrapper for the LAPACK routines zgetrf2 and zgetri.
  subroutine invert_LU_complex_dp(A)
    !> On entry, matrix to invert \( \mathbf{A} \) and
    !> on exit inverted matrix \(\mathbf{A}^{-1}\)
    complex(dp), intent(inout), contiguous :: A(:, :)

    integer, allocatable :: row_pivot_map(:)

    call assert(is_square(A), 'A is not a square matrix.')
    !call assert(abs(determinant_LU(A)) > default_det_tol, 'The absolute value of the determinant of A is zero.')

    allocate(row_pivot_map(size(A, dim=1)))

    call xgetrf2(A, row_pivot_map)
    call xgetri(A, row_pivot_map)
  end subroutine invert_LU_complex_dp


  !> See **[[invert_LU_complex_dp]]**. If this interface is used, the array `A` is unchanged and the matrix
  !> \( \mathbf{A}^{-1} \) is copied to a new array.
  !>
  !> This is a wrapper for the LAPACK routines zgetrf2 and zgetri.
  function inverse_LU_complex_dp(A) result(A_inv)
    !> Matrix to invert \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)

    complex(dp), allocatable :: A_inv(:, :)

    A_inv = A
    call invert_LU_complex_dp(A_inv)
  end function inverse_LU_complex_dp

end module inverse
