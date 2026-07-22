!> Module for hermitian (symmetric) matrix-vector and matrix-matrix multiplication.
!> The interfaces combine wrappers for the LAPACK routines 
!> DSYMV, ZHEMV, DSYMM, ZHEMM,
!> provided as subroutine calls.
module hermitian_matrix_multiplication
  use precision, only: dp, i32
  use constants, only: zone, zzero
#include "asserts.fpp"
  use math_utils, only: is_hermitian
  use lapack_f95_interfaces, only: dsymv, zhemv, dsymm, zhemm

  implicit none

  private
  public :: hermitian_matrix_multiply

  !> Spacing between the elements of an vector.
  integer, parameter :: storage_spacing = 1
  !> Default value for uplo
  character(len=1), parameter :: uplo_default = 'U'
  !> Default value for side
  character(len=1), parameter :: side_default = 'L'
  !> Default value for alpha (double real) for BLAS calls
  real(dp), parameter :: alpha_default_real_dp = 1._dp
  !> Default value for alpha (double complex) for BLAS calls
  complex(dp), parameter :: alpha_default_complex_dp = zone
  !> Default value for beta (double real) for BLAS calls
  real(dp), parameter :: beta_default_real_dp = 0._dp
  !> Default value for beta (double complex) for BLAS calls
  complex(dp), parameter :: beta_default_complex_dp = zzero
  !> Default value for the tolerance used in calls to [[is_hermitian]]
  real(dp), parameter :: tol_default = 1e-10_dp

  !> Calculate the matrix-matrix and matrix-vector product as a subroutine call
  !> \[ \mathbf{C} = \mathbf{A} \cdot \mathbf{B}, \]
  !> \[ \mathbf{c} = \mathbf{A} \cdot \mathbf{b}, \]
  !> where A or B is a hermitian (complex) or symmetric (real) matrix respectively.
  interface hermitian_matrix_multiply
    module procedure :: matrix_vector_multiplication_real_dp, &
                        matrix_vector_multiplication_complex_dp, &
                        matrix_vector_multiplication_real_complex_dp, &
                        matrix_vector_multiplication_complex_real_dp, &
                        matrix_matrix_multiplication_real_dp, &
                        matrix_matrix_multiplication_complex_dp, &
                        matrix_matrix_multiplication_real_complex_dp, &
                        matrix_matrix_multiplication_complex_real_dp
  end interface hermitian_matrix_multiply

contains
  !> Calculate the matrix-vector product between a real symmetric matrix \( \mathbf{A} \)
  !> and a real vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_real_dp(A, b, c, uplo, tol)
    !> Input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    real(dp), intent(out), contiguous :: c(:)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    character(len=1) :: uplo_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a symmetric matrix.')
    CALL_ASSERT(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    CALL_ASSERT(size(A, dim=1) == size(c), 'Number of rows of A needs to be the same as number of elements of c.')

    call dsymv(uplo_, size(A, dim=1), alpha_default_real_dp, A, size(A, dim=1), &
    b, storage_spacing, beta_default_real_dp, C, storage_spacing)
  end subroutine matrix_vector_multiplication_real_dp


  !> Calculate the matrix-vector product between a complex hermitian matrix \( \mathbf{A} \)
  !> and a complex vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_complex_dp(A, b, c, uplo, tol)
    !> Input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    character(len=1) :: uplo_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a hermitian matrix.')
    CALL_ASSERT(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    CALL_ASSERT(size(A, dim=1) == size(c), 'Number of rows of A needs to be the same as number of elements of c.')

    call zhemv(uplo_, size(A, dim=1), alpha_default_complex_dp, A, &
    size(A, dim=1), b, storage_spacing, beta_default_complex_dp, C, storage_spacing)
  end subroutine matrix_vector_multiplication_complex_dp


  !> Calculate the matrix-vector product between a real symmetric matrix \( \mathbf{A} \)
  !> and a complex vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_real_complex_dp(A, b, c, uplo, tol)
    !> Input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    character(len=1) :: uplo_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a hermitian matrix.')
    CALL_ASSERT(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    CALL_ASSERT(size(A, dim=1) == size(c), 'Number of rows of A neecomplexds to be the same as number of elements of c.')

    call zhemv(uplo_, size(A, dim=1), alpha_default_complex_dp, &
    cmplx(A, 0.0_dp, kind=dp), size(A, dim=1), b, storage_spacing, &
    beta_default_complex_dp, C, storage_spacing)
  end subroutine matrix_vector_multiplication_real_complex_dp


  !> Calculate the matrix-vector product between a complex hermitian matrix \( \mathbf{A} \)
  !> and a real vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_complex_real_dp(A, b, c, uplo, tol)
    !> Input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    character(len=1) :: uplo_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a hermitian matrix.')
    CALL_ASSERT(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    CALL_ASSERT(size(A, dim=1) == size(c), 'Number of rows of A needs to be the same as number of elements of c.')

    call zhemv(uplo_, size(A, dim=1), alpha_default_complex_dp, A, &
    size(A, dim=1), cmplx(b, 0.0_dp, kind=dp), storage_spacing, &
    beta_default_complex_dp, C, storage_spacing)
  end subroutine matrix_vector_multiplication_complex_real_dp


! matrix-matrix product

  !> Calculate the matrix-matrix product between two real matrices \( \mathbf{A} \) 
  !> and \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is symmetric.
  subroutine matrix_matrix_multiplication_real_dp(A, B, C, uplo, side, tol)
    !> Input matrices
    real(dp), intent(in), contiguous :: A(:, :), B(:, :)
    !> Output matrix
    real(dp), intent(out), contiguous :: C(:, :)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Define which matrix of \( \mathbf{A} \) and \( \mathbf{B} \) is assumed to be symmetric.
    !>
    !> - for the left matrix: **side** = `'L'` or **side** = `'l'` 
    !>   \( \Rightarrow \mathbf{A} \) is assumed to be symmetric.
    !>
    !> - for the right matrix: **side** = `'R'` or **side** = `'r'` 
    !>   \( \Rightarrow \mathbf{B} \) is assumed to be symmetric.
    !> 
    !> Default is **side** = `'L'`.
    character(len=1), intent(in), optional :: side
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    logical :: is_side_L
    character(len=1) :: uplo_, side_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = side_default
    if (present(side)) side_ = side

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    CALL_ASSERT(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    CALL_ASSERT(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    CALL_ASSERT(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a symmetric matrix.')
      call dsymm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_real_dp, A, size(A, dim=1), B, size(B, dim=1), &
      beta_default_real_dp, C, size(A, dim=1))
    else 
      CALL_ASSERT(is_hermitian(B, tolerance), 'B needs to be a symmetric matrix.')
      call dsymm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_real_dp, B, size(B, dim=2), A, size(A, dim=1), &
      beta_default_real_dp, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_real_dp

  
  !> Calculate the matrix-matrix product between two complex matrices \( \mathbf{A} \)  
  !> and \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is hermitian.
  subroutine matrix_matrix_multiplication_complex_dp(A, B, C, uplo, side, tol)
    !> Input matrices
    complex(dp), intent(in), contiguous :: A(:, :), B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`.
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`.
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Define which matrix of \( \mathbf{A} \) and \( \mathbf{B} \) is assumed to be hermitian.
    !>
    !> - for the left matrix: **side** = `'L'` or **side** = `'l'` 
    !>   \( \Rightarrow \mathbf{A} \) is assumed to be hermitian.
    !>
    !> - for the right matrix: **side** = `'R'` or **side** = `'r'` 
    !>   \( \Rightarrow \mathbf{B} \) is assumed to be hermitian.
    !> 
    !> Default is **side** = `'L'`.
    character(len=1), intent(in), optional :: side
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    logical :: is_side_L
    character(len=1) :: uplo_, side_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = side_default
    if (present(side)) side_ = side

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    CALL_ASSERT(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    CALL_ASSERT(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    CALL_ASSERT(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a hermitian matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_complex_dp, A, size(A, dim=1), B, size(B, dim=1), &
      beta_default_complex_dp, C, size(A, dim=1))
    else
      CALL_ASSERT(is_hermitian(B, tolerance), 'B needs to be a hermitian matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_complex_dp, B, size(B, dim=2), A, size(A, dim=1), &
      beta_default_complex_dp, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_complex_dp


  !> Calculate the matrix-matrix product between a real matrix \( \mathbf{A} \)  
  !> and a complex matrix \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is symmetric or hermitian respectively.
  subroutine matrix_matrix_multiplication_real_complex_dp(A, B, C, uplo, side, tol)
    !> Real input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Complex input matrix \( \mathbf{B} \)
    complex(dp), intent(in), contiguous :: B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`.
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`.
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Define which matrix of \( \mathbf{A} \) and \( \mathbf{B} \) is assumed to be symmetric 
    !> or hermitian respectively.
    !>
    !> - for the left matrix: **side** = `'L'` or **side** = `'l'` 
    !>   \( \Rightarrow \mathbf{A} \) is assumed to be symmetric.
    !>
    !> - for the right matrix: **side** = `'R'` or **side** = `'r'` 
    !>   \( \Rightarrow \mathbf{B} \) is assumed to be hermitian.
    !> 
    !> Default is **side** = `'L'`.
    character(len=1), intent(in), optional :: side
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    logical :: is_side_L
    character(len=1) :: uplo_, side_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = side_default
    if (present(side)) side_ = side

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    CALL_ASSERT(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    CALL_ASSERT(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    CALL_ASSERT(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a symmetric matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_complex_dp, cmplx(A, 0.0_dp, kind=dp), size(A, dim=1), &
      B, size(B, dim=1), beta_default_complex_dp, C, size(A, dim=1))
    else 
      CALL_ASSERT(is_hermitian(B, tolerance), 'B needs to be a hermitian matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_complex_dp, B, size(B, dim=2), cmplx(A, 0.0_dp, kind=dp), &
      size(A, dim=1), beta_default_complex_dp, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_real_complex_dp


  !> Calculate the matrix-matrix product between a complex matrix \( \mathbf{A} \)  
  !> and a real matrix \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is hermitian or symmetric respectively.
  subroutine matrix_matrix_multiplication_complex_real_dp(A, B, C, uplo, side, tol)
    !> Real input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Complex input matrix \( \mathbf{B} \)
    real(dp), intent(in), contiguous :: B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Define if the upper or the lower triangular part of 
    !> the hermitian matrix \( \mathbf{A} \) is to be referenced:
    !>
    !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`.
    !> 
    !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`.
    !>
    !> Default is **uplo** = `'U'`.
    character(len=1), intent(in), optional :: uplo
    !> Define which matrix of \( \mathbf{A} \) and \( \mathbf{B} \) is assumed to be symmetric 
    !> or hermitian respectively.
    !>
    !> - for the left matrix: **side** = `'L'` or **side** = `'l'` 
    !>   \( \Rightarrow \mathbf{A} \) is assumed to be hermitian.
    !>
    !> - for the right matrix: **side** = `'R'` or **side** = `'r'` 
    !>   \( \Rightarrow \mathbf{B} \) is assumed to be symmetric.
    !> 
    !> Default is **side** = `'L'`.
    character(len=1), intent(in), optional :: side
    !> Tolerance for checking if the matrix is hermitian
    real(dp), intent(in), optional :: tol

    logical :: is_side_L
    character(len=1) :: uplo_, side_
    real(dp) ::  tolerance

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = side_default
    if (present(side)) side_ = side

    tolerance = tol_default
    if ( present(tol) ) tolerance = tol

    CALL_ASSERT(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    CALL_ASSERT(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    CALL_ASSERT(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    CALL_ASSERT(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    CALL_ASSERT(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      CALL_ASSERT(is_hermitian(A, tolerance), 'A needs to be a symmetric matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_complex_dp, A, size(A, dim=1), cmplx(B, 0.0_dp, kind=dp), &
      size(B, dim=1), beta_default_complex_dp, C, size(A, dim=1))
    else 
      CALL_ASSERT(is_hermitian(B, tolerance), 'B needs to be a symmetric matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), &
      alpha_default_complex_dp, cmplx(B, 0.0_dp, kind=dp), size(B, dim=2), A, &
      size(A, dim=1), beta_default_complex_dp, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_complex_real_dp

end module hermitian_matrix_multiplication