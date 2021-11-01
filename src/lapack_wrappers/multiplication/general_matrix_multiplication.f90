!> Module for general matrix-vector and matrix-matrix multiplication.
!> The interfaces combine wrappers for the LAPACK routines 
!> DGEMV, ZGEMVM, DGEMM, ZGEMM.
module general_matrix_multiplication
  use precision, only: dp
  use constants, only: zone, zzero
  use lapack_f95_interfaces, only: dgemv, zgemv, dgemm, zgemm
  use asserts, only: assert

  implicit none
  
  private
  public :: matrix_multiply

  !> Calculate the matrix-vector, vector-matrix and matrix-matrix, product for general matrices 
  !> with a subroutine call:
  !> \[ \mathbf{c} = \mathbf{A} \cdot \mathbf{b}, \]
  !> \[ \mathbf{c} = \mathbf{a} \cdot \mathbf{B}, \]
  !> \[ \mathbf{C} = \mathbf{A} \cdot \mathbf{B}. \]
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
  !> where op takes into account that the matrix can be used also as transposed (see **trans_A**). 
  subroutine matrix_vector_multiplication_real_dp(A, b, c, trans_A) 
    !> Input matrix     
    real(dp), intent(in), contiguous :: A(:,:)
    !> Input vector  
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    real(dp), intent(out), contiguous :: c(:)
    !> Define if \( \mathbf{A} \) is used as is for 'N', 'n' or as transposed for 'T', 't'.
    !> Default is 'N'. 
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = 'N'
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call dgemv(trans_A_, size(A, dim=1), size(A, dim=2), 1.0_dp, A, size(A, dim=1), b, 1, 0.0_dp, c, 1)
  end subroutine matrix_vector_multiplication_real_dp


  !> Calculate the matrix-vector product between a complex matrix \(\mathbf{A}\) and a complex vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where op takes into account that the matrix can be used also as transposed or conjugate-transposed (see **trans_A**). 
  subroutine matrix_vector_multiplication_complex_dp(A, b, c, trans_A)
    !> Input matrix 
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector  
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define if \( \mathbf{A} \) is used as is for 'N', 'n', as transposed for 'T', 't' or 
    !> as complex conjugate-transpose for 'c', 'c'.
    !> Default is 'N'.
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = 'N'
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call zgemv(trans_A_, size(A, dim=1), size(A, dim=2), zone, A, size(A, dim=1), b, 1, zzero, c, 1)
  end subroutine matrix_vector_multiplication_complex_dp


  !> Calculate the matrix-vector product between a real matrix \(\mathbf{A}\) and a complex vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where op takes into account that the matrix can be used also as transposed or conjugate-transposed (see **trans_A**). 
  subroutine matrix_vector_multiplication_real_complex_dp(A, b, c, trans_A)
    !> Input matrix 
    real(dp), intent(in), contiguous :: A(:, :)
    !> Input vector  
    complex(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define if \( \mathbf{A} \) is used as is for 'N', 'n', as transposed for 'T', 't' or 
    !> as complex conjugate-transpose for 'c', 'c'.
    !> Default is 'N'.
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = 'N'
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call zgemv(trans_A_, size(A, dim=1), size(A, dim=2), zone, zone * A, size(A, dim=1), b, 1, zzero, c, 1)
  end subroutine matrix_vector_multiplication_real_complex_dp


  !> Calculate the matrix-vector product between a complex matrix \(\mathbf{A}\) and a real vector \(\mathbf{b}\):
  !> \[
  !>     \mathbf{c} = \text{op}(\mathbf{A}) \cdot \mathbf{b},
  !> \]
  !> where op takes into account that the matrix can be used also as transposed or conjugate-transposed (see **trans_A**). 
  subroutine matrix_vector_multiplication_complex_real_dp(A, b, c, trans_A)
    !> Input matrix 
    complex(dp), intent(in), contiguous :: A(:, :)
    !> Input vector  
    real(dp), intent(in), contiguous :: b(:)
    !> Output vector
    complex(dp), intent(out), contiguous :: c(:)
    !> Define if \( \mathbf{A} \) is used as is for 'N', 'n', as transposed for 'T', 't' or 
    !> as complex conjugate-transpose for 'c', 'c'.
    !> Default is 'N'.
    character(len=1), intent(in), optional :: trans_A                    
    
    character(len=1) :: trans_A_

    trans_A_ = 'N'
    if (present(trans_A)) trans_A_ = trans_A
    
    call assert(size(c) == integer_gemv(shape(A), size(b), trans_A_), &
               'The number of elements of c must be equal the number of rows of op(A).')

    call zgemv(trans_A_, size(A, dim=1), size(A, dim=2), zone, A, size(A, dim=1), zone * b, 1, zzero, c, 1)
  end subroutine matrix_vector_multiplication_complex_real_dp


! matrix-matrix product

  !> Calculate the matrix-matrix product between two real matrices \(\mathbf{A}\) and \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where op takes into account that the matrices can be used also as transposed (see **trans_A** and **trans_B**).
  subroutine matrix_matrix_multiplication_real_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    real(dp), intent(in), contiguous :: A(:, :), B(:, :)
    !> Output matrix 
    real(dp), intent(out), contiguous :: C(:, :)
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n')
    !> or the transposed is used ('T' or 't').
    character(len=1), intent(in), optional :: trans_A, trans_B

    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = 'N'
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = 'N'
    if(present(trans_B)) trans_B_ = trans_B
               
    call integer_gemm(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call dgemm(trans_A_, trans_B_, M, N, K, 1.0_dp, A, LDA, B, LDB, 0.0_dp, C, LDC)
  end subroutine matrix_matrix_multiplication_real_dp


  !> Calculate the matrix-matrix product between two complex matrices \(\mathbf{A}\) and \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where op takes into account that the matrices can be also used as transposed or 
  !> conjugate-transposed (see **trans_A** and **trans_B**). 
  subroutine matrix_matrix_multiplication_complex_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    complex(dp), intent(in), contiguous  :: A(:, :), B(:, :)
    !> Output matrix
    complex(dp), intent(out), contiguous  :: C(:, :)
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n')
    !> or the transposed is used ('T' or 't').
    character(len=1), intent(in), optional :: trans_A, trans_B
    
    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = 'N'
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = 'N'
    if(present(trans_B)) trans_B_ = trans_B
               
    call integer_gemm(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call zgemm(trans_A_, trans_B_, M, N, K, zone, A, LDA, B, LDB, zzero, C, LDC)
  end subroutine matrix_matrix_multiplication_complex_dp

  
  !> Calculate the matrix-matrix product between a real matrix \(\mathbf{A}\) and a real matrix \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where op takes into account that the matrices can be also used as transposed or 
  !> conjugate-transposed (see **trans_A** and **trans_B**). 
  subroutine matrix_matrix_multiplication_real_complex_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    real(dp), intent(in), contiguous :: A(:, :)
    complex(dp), intent(in), contiguous :: B(:, :)
    !> Output matrix 
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n')
    !> or the transposed is used ('T' or 't').
    character(len=1), intent(in), optional :: trans_A, trans_B

    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = 'N'
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = 'N'
    if(present(trans_B)) trans_B_ = trans_B
               
    call integer_gemm(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call zgemm(trans_A_, trans_B_, M, N, K, zone, zone * A, LDA, B, LDB, zzero, C, LDC)
  end subroutine matrix_matrix_multiplication_real_complex_dp


  !> Calculate the matrix-matrix product between a complex matrix \(\mathbf{A}\) and a real matrix \(\mathbf{B}\).
  !> \[
  !>    \mathbf{C} = \text{op}(\mathbf{A}) \cdot \text{op}(\mathbf{B}),
  !> \]
  !> where op takes into account that the matrices can be also used as transposed or 
  !> conjugate-transposed (see **trans_A** and **trans_B**). 
  subroutine matrix_matrix_multiplication_complex_real_dp(A, B, C, trans_A, trans_B)
    !> Input matrices 
    complex(dp), intent(in), contiguous :: A(:, :)
    real(dp), intent(in), contiguous :: B(:, :)
    !> Output matrix 
    complex(dp), intent(out), contiguous :: C(:, :)
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n')
    !> or the transposed is used ('T' or 't').
    character(len=1), intent(in), optional :: trans_A, trans_B

    integer :: M, N, K, LDA, LDB, LDC
    character(len=1) :: trans_A_, trans_B_
    
    trans_A_ = 'N'
    if(present(trans_A)) trans_A_ = trans_A

    trans_B_ = 'N'
    if(present(trans_B)) trans_B_ = trans_B
               
    call integer_gemm(shape(A), shape(B), trans_A_, trans_B_,  M, N, K, LDA, LDB, LDC)

    call assert(all(shape(C) == [M, N]), 'The number of rows of C must be equal to the number of rows of op(A) &
               and the number of culomns must be equal to the number of columns of op(B).')

    call zgemm(trans_A_, trans_B_, M, N, K, zone, A, LDA, zone * B, LDB, zzero, C, LDC)
  end subroutine matrix_matrix_multiplication_complex_real_dp
  
  
  ! Setup functions for LAPACK routine call

  !> Setup the integers for *gemv calls.
  integer function integer_gemv(shape_A, size_B, trans_A)
    !> Shape of the input arrays
    integer, intent(in) :: shape_A(2), size_B
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n'),
    !> the transposed ('T' or 't') or the complex conjugated ('C', 'c') is used.
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
  subroutine integer_gemm(shape_A, shape_B, trans_A, trans_B, M, N, K, LDA, LDB, LDC)
    !> Shape of the input arrays
    integer, intent(in) :: shape_A(2), shape_B(2)
    !> Determine if \( \mathbf{A} \) or \( \mathbf{B} \) are used as is ('N' or 'n'),
    !> the transposed ('T' or 't') or the complex conjugated ('C', 'c') is used.
    character(len=1), intent(in) :: trans_A, trans_B
    !> Integers for *func_matrix_multiply call
    integer, intent(out) :: M, N, K, LDA, LDB, LDC

    logical :: trans_A_, trans_B_

    trans_A_ = any(trans_A == ["T", "t", "C", "c"])
    trans_B_ = any(trans_B == ["T", "t", "C", "c"])

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
  end subroutine integer_gemm


end module general_matrix_multiplication