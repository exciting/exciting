!> Module for general matrix-vector and matrix-matrix multiplication.
!> The interfaces combine wrappers for the LAPACK routines 
!> DGEMV, ZGEMV, DGEMM, ZGEMM.
module general_matrix_multiplication
  use precision, only: dp
  use constants, only: zone, zzero
  use lapack_f95_interfaces, only: dgemv, zgemv, dgemm, zgemm
  use asserts, only: assert

  implicit none
  
  private
  public :: matrix_multiply, gemm_parameters


  !> Default for **trans_A** and **trans_B** respectively which define if the matrix is used, its transpose or its
  !> conjugate transpose in case of complex matrices. Default is to use the matrix (`'N'`).
  character(len=1), parameter :: default_trans_char = 'N'


  !> Calculate for general matrices the matrix-vector and the matrix-matrix product.
  !>
  !> Capital letters denote matrices and small letters vectors.
  !> The matrix-vector product
  !> \[ \mathbf{c} = \mathbf{A} \cdot \mathbf{b} \]
  !> is computed with the subroutine call `call matrix_multiply(A, b, c)`.
  !>
  !> The matrix-matrix product
  !> \[ \mathbf{C} = \mathbf{A} \cdot \mathbf{B} \]
  !> is computed with the subroutine call `call matrix_multiply(A, B, C)`.
  interface matrix_multiply
    module procedure :: matrix_vector_multiplication_real_dp, &
                        matrix_vector_multiplication_complex_dp, &
                        matrix_vector_multiplication_real_complex_dp, &
                        matrix_vector_multiplication_complex_real_dp, &
                        matrix_matrix_multiplication_real_dp,&
                        matrix_matrix_multiplication_complex_dp, &
                        matrix_matrix_multiplication_real_complex_dp, &
                        matrix_matrix_multiplication_complex_real_dp
  end interface matrix_multiply

contains


! matrix-vector product

  !> Calculate the matrix-vector product between a real matrix \(\mathbf{A}\) and a real vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`.
  subroutine matrix_vector_multiplication_real_dp(A, b, c, trans_A) 
    !> Input matrix     
    real(dp), intent(in), contiguous :: A(:,:)
    !> Input vector  
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    real(dp), intent(out), contiguous :: c(:)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = default_trans_char
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call dgemv(trans_A_, size(A, dim=1), size(A, dim=2), 1.0_dp, A, size(A, dim=1), b, 1, 0.0_dp, c, 1)
  end subroutine matrix_vector_multiplication_real_dp


  !> Calculate the matrix-vector product between a complex matrix \(\mathbf{A}\) and a complex vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \) if **trans_A** = `'C'` or **trans_A** = `'c'`.
  subroutine matrix_vector_multiplication_complex_dp(A, b, c, trans_A)
    !> Input matrix 
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector  
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> - if **trans_A** = `'C'` or **trans_A** = `'c'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = default_trans_char
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call zgemv(trans_A_, size(A, dim=1), size(A, dim=2), zone, A, size(A, dim=1), b, 1, zzero, c, 1)
  end subroutine matrix_vector_multiplication_complex_dp


  !> Calculate the matrix-vector product between a real matrix \(\mathbf{A}\) and a complex vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`.
  subroutine matrix_vector_multiplication_real_complex_dp(A, b, c, trans_A)
    !> Input matrix 
    real(dp), intent(in), contiguous :: A(:, :)
    !> Input vector  
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = default_trans_char
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call zgemv(trans_A_, size(A, dim=1), size(A, dim=2), zone, cmplx(A, 0.0_dp, kind=dp), size(A, dim=1), b, 1, zzero, c, 1)
  end subroutine matrix_vector_multiplication_real_complex_dp


  !> Calculate the matrix-vector product between a complex matrix \(\mathbf{A}\) and a real vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \) if **trans_A** = `'C'` or **trans_A** = `'c'`.
  subroutine matrix_vector_multiplication_complex_real_dp(A, b, c, trans_A)
    !> Input matrix 
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector  
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> - if **trans_A** = `'C'` or **trans_A** = `'c'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = default_trans_char
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call zgemv(trans_A_, size(A, dim=1), size(A, dim=2), zone, A, size(A, dim=1), cmplx(b, 0.0_dp, kind=dp), 1, zzero, c, 1)
  end subroutine matrix_vector_multiplication_complex_real_dp


! matrix-matrix product

  !> Calculate the matrix-matrix product between two real matrices \(\mathbf{A}\) and \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`.
  !>
  !> and \( \text{op}(\mathbf{B}) \) is, depending on **trans_B**, one of 
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B} \) if **trans_B** = `'N'` or **trans_B** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \) if **trans_B** = `'T'` or **trans_B** = `'t'`.
  subroutine matrix_matrix_multiplication_real_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    real(dp), intent(in), contiguous :: A(:, :), B(:, :)
    !> Output matrix 
    real(dp), intent(out), contiguous :: C(:, :)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A
    !> Define \(  \text{op}(\mathbf{B}) \): 
    !> 
    !> - if **trans_B** = `'N'` or **trans_B** = `'n'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B} \).
    !>
    !> - if **trans_B** = `'T'` or **trans_B** = `'t'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_B

    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = default_trans_char
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = default_trans_char
    if(present(trans_B)) trans_B_ = trans_B
               
    call gemm_parameters(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call dgemm(trans_A_, trans_B_, M, N, K, 1.0_dp, A, LDA, B, LDB, 0.0_dp, C, LDC)
  end subroutine matrix_matrix_multiplication_real_dp


  !> Calculate the matrix-matrix product between two complex matrices \(\mathbf{A}\) and \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \) if **trans_A** = `'C'` or **trans_A** = `'c'`.
  !>
  !> and \( \text{op}(\mathbf{B}) \) is, depending on **trans_B**, one of 
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B} \) if **trans_B** = `'N'` or **trans_B** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \) if **trans_B** = `'T'` or **trans_B** = `'t'`,
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B}^\dagger \) if **trans_B** = `'C'` or **trans_B** = `'c'`.
  subroutine matrix_matrix_multiplication_complex_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    complex(dp), intent(in), contiguous  :: A(:, :), B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous  :: C(:, :)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> - if **trans_A** = `'C'` or **trans_A** = `'c'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A
    !> Define \(  \text{op}(\mathbf{B}) \): 
    !> 
    !> - if **trans_B** = `'N'` or **trans_B** = `'n'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B} \).
    !>
    !> - if **trans_B** = `'T'` or **trans_B** = `'t'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \).
    !>
    !> - if **trans_B** = `'C'` or **trans_B** = `'c'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\dagger \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_B
    
    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = default_trans_char
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = default_trans_char
    if(present(trans_B)) trans_B_ = trans_B
               
    call gemm_parameters(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call zgemm(trans_A_, trans_B_, M, N, K, zone, A, LDA, B, LDB, zzero, C, LDC)
  end subroutine matrix_matrix_multiplication_complex_dp

  
  !> Calculate the matrix-matrix product between a real matrix \(\mathbf{A}\) and a real matrix \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`.
  !>
  !> and \( \text{op}(\mathbf{B}) \) is, depending on **trans_B**, one of 
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B} \) if **trans_B** = `'N'` or **trans_B** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \) if **trans_B** = `'T'` or **trans_B** = `'t'`,
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B}^\dagger \) if **trans_B** = `'C'` or **trans_B** = `'c'`.
  subroutine matrix_matrix_multiplication_real_complex_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    real(dp), intent(in), contiguous :: A(:, :)
    complex(dp), intent(in), contiguous :: B(:, :)
    !> Output matrix 
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A
    !> Define \(  \text{op}(\mathbf{B}) \): 
    !> 
    !> - if **trans_B** = `'N'` or **trans_B** = `'n'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B} \).
    !>
    !> - if **trans_B** = `'T'` or **trans_B** = `'t'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \).
    !>
    !> - if **trans_B** = `'C'` or **trans_B** = `'c'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\dagger \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_B

    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = default_trans_char
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = default_trans_char
    if(present(trans_B)) trans_B_ = trans_B
               
    call gemm_parameters(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call zgemm(trans_A_, trans_B_, M, N, K, zone, cmplx(A, 0.0_dp, kind=dp), LDA, B, LDB, zzero, C, LDC)
  end subroutine matrix_matrix_multiplication_real_complex_dp


  !> Calculate the matrix-matrix product between a complex matrix \(\mathbf{A}\) and a real matrix \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where \( \text{op}(\mathbf{A}) \) is, depending on **trans_A**, one of 
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A} \) if **trans_A** = `'N'` or **trans_A** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\top \) if **trans_A** = `'T'` or **trans_A** = `'t'`,
  !>
  !> - \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \) if **trans_A** = `'C'` or **trans_A** = `'c'`.
  !>
  !> and \( \text{op}(\mathbf{B}) \) is, depending on **trans_B**, one of 
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B} \) if **trans_B** = `'N'` or **trans_B** = `'n'`,
  !>
  !> - \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \) if **trans_B** = `'T'` or **trans_B** = `'t'`.
  subroutine matrix_matrix_multiplication_complex_real_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    complex(dp), intent(in), contiguous :: A(:, :)
    real(dp), intent(in), contiguous :: B(:, :)
    !> Output matrix 
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Define \(  \text{op}(\mathbf{A}) \): 
    !> 
    !> - if **trans_A** = `'N'` or **trans_A** = `'n'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A} \).
    !>
    !> - if **trans_A** = `'T'` or **trans_A** = `'t'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \).
    !>
    !> - if **trans_A** = `'C'` or **trans_A** = `'c'`,
    !>   \( \text{op}(\mathbf{A}) = \mathbf{A}^\dagger \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_A
    !> Define \(  \text{op}(\mathbf{B}) \): 
    !> 
    !> - if **trans_B** = `'N'` or **trans_B** = `'n'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B} \).
    !>
    !> - if **trans_B** = `'T'` or **trans_B** = `'t'`,
    !>   \( \text{op}(\mathbf{B}) = \mathbf{B}^\top \).
    !>
    !> Default is gven by [[default_trans_char]].
    character(len=1), intent(in), optional :: trans_B


    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = default_trans_char
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = default_trans_char
    if(present(trans_B)) trans_B_ = trans_B
               
    call gemm_parameters(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call zgemm(trans_A_, trans_B_, M, N, K, zone, A, LDA, cmplx(B, 0.0_dp, kind=dp), LDB, zzero, C, LDC)
  end subroutine matrix_matrix_multiplication_complex_real_dp


  ! Setup functions for LAPACK routine call

  !> Setup the integers for *gemv calls.
  integer function integer_gemv(shape_A, size_B, trans_A)
    !> Shape of the input arrays
    integer, intent(in) :: shape_A(2), size_B
    !> Determine if \( \mathbf{A} \) is used as is ('N' or 'n'),
    !> the transpose (`'T'` or `'t'`) or the conjugate transpose (`'C'`, `'c'`) is used.
    character(len=1), intent(in) :: trans_A
    !> Integers for *func_matrix_multiply call

    if (any(trans_A == ['N', 'n'])) then
      call assert(shape_A(2) == size_B, &
                 'Number of columns of A must be the same as number of elements of B.')
      integer_gemv = shape_A(1)

    else if (any(trans_A == ['T', 't', 'C', 'c'])) then
      call assert(shape_A(1) == size_B, &
                 'Number of columns of A**T must be the same as number of elements of B.')
      integer_gemv = shape_A(2)

    else
      call assert(.false., &
                 'trans_A must be one of "N", "n", "T", "t", "C", "c".')
    end if
 
  end function integer_gemv

  !> Setup the integers for *gemm calls.
  subroutine gemm_parameters(shape_A, shape_B, trans_A, trans_B, M, N, K, LDA, LDB, LDC)
    !> Shape of the input arrays
    integer, intent(in) :: shape_A(2), shape_B(2)
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n'),
    !> the transposed (`'T'` or `'t'`) or the complex conjugated (`'C'`, `'c'`) is used.
    character(len=1), intent(in) :: trans_A, trans_B
    !> Integers for *gemm calls
    integer, intent(out) :: M, N, K, LDA, LDB, LDC

    logical :: trans_A_, trans_B_

    trans_A_ = any(trans_A == ['T', 't', 'C', 'c'])
    trans_B_ = any(trans_B == ['T', 't', 'C', 'c'])

    if ((.not. trans_A_) .and. (.not. trans_B_)) then
      call assert(shape_A(2) == shape_B(1), &
                  'Number of columns of A must be the same as number of rows of B.')

      M = shape_A(1); K = shape_A(2); N = shape_B(2)
      LDA = M; LDB = K; LDC = M

    else if ((.not. trans_A_) .and. trans_B_) then
      call assert(shape_A(2) == shape_B(2), &
                  'Number of columns of A must be the same as number of rows of B**T.')

      M = shape_A(1); K = shape_A(2); N = shape_B(1)
      LDA = M; LDB = N; LDC = M

    else if (trans_A_ .and. (.not. trans_B_)) then
      call assert(shape_A(1) == shape_B(1), &
                 'Number of columns of A**T must be the same as number of rows of B.')

      M = shape_A(2); K = shape_A(1); N = shape_B(2)
      LDA = K; LDB = K; LDC = M

    else if (trans_A_  .and. trans_B_) then
      call assert(shape_A(1) == shape_B(2), &
                 'Number of columns of A**T must be the same as number of rows of B**T.')

      M = shape_A(2);  K = shape_A(1); N = shape_B(1)
      LDA = K; LDB = N; LDC = M

    else
      call assert(.false., 'trans_A must be one of "N", "n", "T", "t" "C", "c".')
    end if
  end subroutine gemm_parameters

end module general_matrix_multiplication