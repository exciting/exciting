!> Module for hermitian (symmetric) matrix-vector and matrix-matrix multiplication.
!> The interfaces combine wrappers for the LAPACK routines 
!> DSYMV, ZHEMV, DSYMM, ZHEMM,
!> provided as subroutine calls and function calls.
!> Use the subroutine call whenever efficiency is a major issue or 
!> you want to mutate a matrix or a slice of an array.
module hermitian_matrix_multiplication
  use precision, only: dp
  use constants, only: zone, zzero
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use math_utils, only: is_hermitian
  use lapack_f95_interfaces, only: dsymv, zhemv, dsymm, zhemm

  implicit none

  private
  public :: hermitian_matrix_multiply

  !> Spacing between the elements of an vector.
  integer, parameter :: storage_spacing = 1
  !> Default value for uplo
  character(len=1), parameter :: uplo_default = 'U'


  !> Calculate the matrix-matrix, matrix-vector and vector-matrix product as a subroutine call
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
  subroutine matrix_vector_multiplication_real_dp(A, b, c, uplo)
    !> Input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    real(dp), intent(out), contiguous :: c(:)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo

    character(len=1) :: uplo_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(is_hermitian(A), 'A needs to be a symmetric matrix.')
    call assert(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    call assert(size(A, dim=1) == size(c), 'Number of rows of A needs to be the same as number of elements of c.')

    call dsymv(uplo_, size(A, dim=1), 1.0_dp, A, size(A, dim=1), b, storage_spacing, 0.0_dp, C, storage_spacing)
  end subroutine matrix_vector_multiplication_real_dp


  !> Calculate the matrix-vector product between a complex hermitian matrix \( \mathbf{A} \)
  !> and a complex vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_complex_dp(A, b, c, uplo)
    !> Input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo

    character(len=1) :: uplo_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(is_hermitian(A), 'A needs to be a hermitian matrix.')
    call assert(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    call assert(size(A, dim=1) == size(c), 'Number of rows of A needs to be the same as number of elements of c.')
    
    call zhemv(uplo_, size(A, dim=1), zone, A, size(A, dim=1), b, storage_spacing, zzero, C, storage_spacing)
  end subroutine matrix_vector_multiplication_complex_dp


  !> Calculate the matrix-vector product between a real symmetric matrix \( \mathbf{A} \)
  !> and a complex vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_real_complex_dp(A, b, c, uplo)
    !> Input matrix \( \mathbf{A} \)
    real(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo

    character(len=1) :: uplo_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(is_hermitian(A), 'A needs to be a hermitian matrix.')
    call assert(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    call assert(size(A, dim=1) == size(c), 'Number of rows of A neecomplexds to be the same as number of elements of c.')
    
    call zhemv(uplo_, size(A, dim=1), zone, cmplx(A, 0.0_dp, kind=dp), size(A, dim=1), b, storage_spacing, zzero, C, storage_spacing)
  end subroutine matrix_vector_multiplication_real_complex_dp


  !> Calculate the matrix-vector product between a complex hermitian matrix \( \mathbf{A} \)
  !> and a real vector \( \mathbf{b} \):
  !> \[
  !>    \mathbf{c} = \mathbf{A} \cdot \mathbf{b}.
  !> \]
  subroutine matrix_vector_multiplication_complex_real_dp(A, b, c, uplo)
    !> Input matrix \( \mathbf{A} \)
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector \( \mathbf{b} \)
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo

    character(len=1) :: uplo_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo 

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(is_hermitian(A), 'A needs to be a hermitian matrix.')
    call assert(size(A, dim=2) == size(b), 'Number of columns of A needs to be the same as number of elements of b.')
    call assert(size(A, dim=1) == size(c), 'Number of rows of A needs to be the same as number of elements of c.')
    
    call zhemv(uplo_, size(A, dim=1), zone, A, size(A, dim=1), cmplx(b, 0.0_dp, kind=dp), storage_spacing, zzero, C, storage_spacing)
  end subroutine matrix_vector_multiplication_complex_real_dp


! matrix-matrix product

  !> Calculate the matrix-matrix product between two real matrices \( mathbf{A} \) 
  !> and \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is symmetric.
  subroutine matrix_matrix_multiplication_real_dp(A, B, C, uplo, side)
    !> Input matrices
    real(dp), intent(in), contiguous :: A(:, :), B(:, :)
    !> Output matrix
    real(dp), intent(out), contiguous :: C(:, :)
    !> Specify whether the upper or lower triangular part of 
    !> the symmetric matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo
    !> Specify from which which of \( mathbf{A} \) and \( \mathbf{B} \) is symmetric.
    !> 'L' or 'l' for left:  \( mathbf{A} \) is symmetric (default)
    !> 'R' or 'r' for right: \( mathbf{B} \) is symmetric
    character(len=1), intent(in), optional :: side

    logical :: is_side_L
    character(len=1) :: uplo_, side_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = 'L'
    if (present(side)) side_ = side

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    call assert(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    call assert(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    call assert(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      call assert(is_hermitian(A), 'A needs to be a symmetric matrix.')
      call dsymm(side_, uplo_, size(A, dim=1), size(B, dim=2), 1.0_dp, A, size(A, dim=1), B, size(B, dim=1), 0.0_dp, C, size(A, dim=1))
    else 
      call assert(is_hermitian(B), 'B needs to be a symmetric matrix.')
      call dsymm(side_, uplo_, size(A, dim=1), size(B, dim=2), 1.0_dp, B, size(B, dim=2), A, size(A, dim=1), 0.0_dp, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_real_dp

  
  !> Calculate the matrix-matrix product between two complex matrices \( mathbf{A} \)  
  !> and \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is hermitian.
  subroutine matrix_matrix_multiplication_complex_dp(A, B, C, uplo, side)
    !> Input matrices
    complex(dp), intent(in), contiguous :: A(:, :), B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo
    !> Specify from which which of \( mathbf{A} \) and \( \mathbf{B} \) is hermitian.
    !> 'L' or 'l' for left:  \( mathbf{A} \) is hermitian (default)
    !> 'R' or 'r' for right: \( mathbf{B} \) is hermitian
    character(len=1), intent(in), optional :: side

    logical :: is_side_L
    character(len=1) :: uplo_, side_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = 'L'
    if (present(side)) side_ = side

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    call assert(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    call assert(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    call assert(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      call assert(is_hermitian(A), 'A needs to be a hermitian matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), zone, A, size(A, dim=1), B, size(B, dim=1), zzero, C, size(A, dim=1))
    else
      call assert(is_hermitian(B), 'B needs to be a hermitian matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), zone, B, size(B, dim=2), A, size(A, dim=1), zzero, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_complex_dp


  !> Calculate the matrix-matrix product between a real matrix \( mathbf{A} \)  
  !> and a complex matrix \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is symmetric or hermitian respectively.
  subroutine matrix_matrix_multiplication_real_complex_dp(A, B, C, uplo, side)
    !> Real input matrix \( \mathbf{A} \)
    real(dp), intent(in) :: A(:, :)
    !> Complex input matrix \( \mathbf{B} \)
    complex(dp), intent(in) :: B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo
    !> Specify from which which of \( mathbf{A} \) and \( \mathbf{B} \) is hermitian.
    !> 'L' or 'l' for left:  \( mathbf{A} \) is hermitian (default)
    !> 'R' or 'r' for right: \( mathbf{B} \) is hermitian
    character(len=1), intent(in), optional :: side

    logical :: is_side_L
    character(len=1) :: uplo_, side_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = 'L'
    if (present(side)) side_ = side

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    call assert(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    call assert(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    call assert(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      call assert(is_hermitian(A), 'A needs to be a symmetric matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), zone, cmplx(A, 0.0_dp, kind=dp), size(A, dim=1), B, size(B, dim=1), zzero, C, size(A, dim=1))
    else 
      call assert(is_hermitian(B), 'B needs to be a hermitian matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), zone, B, size(B, dim=2), cmplx(A, 0.0_dp, kind=dp), size(A, dim=1), zzero, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_real_complex_dp

  !> Calculate the matrix-matrix product between a complex matrix \( mathbf{A} \)  
  !> and a real matrix \( \mathbf{B} \)
  !> \[
  !>    \mathbf{C} = \mathbf{A} \cdot \mathbf{B},
  !> \]
  !> where one of both is hermitian or symmetric respectively.
  subroutine matrix_matrix_multiplication_complex_real_dp(A, B, C, uplo, side)
    !> Real input matrix \( \mathbf{A} \)
    complex(dp), intent(in) :: A(:, :)
    !> Complex input matrix \( \mathbf{B} \)
    real(dp), intent(in) :: B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Specify whether the upper or lower triangular part of 
    !> the hermitian matrix A is to be referenced
    !> 'U' or 'u' upper triangular part (default)
    !> 'L' or 'l' lower triangular part
    character(len=1), intent(in), optional :: uplo
    !> Specify from which which of \( mathbf{A} \) and \( \mathbf{B} \) is hermitian.
    !> 'L' or 'l' for left:  \( mathbf{A} \) is hermitian (default)
    !> 'R' or 'r' for right: \( mathbf{B} \) is hermitian
    character(len=1), intent(in), optional :: side

    logical :: is_side_L
    character(len=1) :: uplo_, side_

    uplo_ = uplo_default
    if (present(uplo)) uplo_ = uplo

    side_ = 'L'
    if (present(side)) side_ = side

    call assert(any(uplo_ == ['U', 'u', 'L', 'l']), 'uplo needs to be one of "U", "u", "L", "l".')
    call assert(any(side_ == ['L', 'l', 'R', 'r']), 'side needs to be one of "L", "l", "R" or "r".')
    call assert(size(A, dim=2) == size(B, dim=1), 'Number of columns of A needs to be the same as the number of rows of B.')
    call assert(size(C, dim=1) == size(A, dim=1), 'The number of rows of C must be equal to the number of rows of A.')
    call assert(size(C, dim=2) == size(B, dim=2), 'The number of columns of C must be equal to the number of columns of B.')

    is_side_L = any(side_ == ['L', 'l'])

    if (is_side_L) then
      call assert(is_hermitian(A), 'A needs to be a symmetric matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), zone, A, size(A, dim=1), cmplx(B, 0.0_dp, kind=dp), size(B, dim=1), zzero, C, size(A, dim=1))
    else 
      call assert(is_hermitian(B), 'B needs to be a symmetric matrix.')
      call zhemm(side_, uplo_, size(A, dim=1), size(B, dim=2), zone, cmplx(B, 0.0_dp, kind=dp), size(B, dim=2), A, size(A, dim=1), zzero, C, size(A, dim=1))
    end if
  end subroutine matrix_matrix_multiplication_complex_real_dp

end module hermitian_matrix_multiplication