!> Unit tests for general matrix multiplication.
module general_matrix_multiplication_test
  use precision, only: dp
  use constants, only: zi
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close
  use mock_arrays, only: real_vector_5, &
                         complex_vector_5, &
                         real_vector_7, &
                         complex_vector_7, &
                         real_matrix_5x7, &
                         real_matrix_7x5, &
                         complex_matrix_5x7, &
                         complex_matrix_7x5
  use general_matrix_multiplication, only: matrix_multiply

  implicit none

  private
  public :: general_matrix_multiplication_test_driver
  

contains


  !> Run tests for general matrix multiplication
  subroutine general_matrix_multiplication_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 35

    call test_report%init(n_assertions, mpiglobal)

    ! Run unit tests

    call test_real_matrix_vector_multiplication(test_report)

    call test_complex_matrix_vector_multiplication(test_report)

    call test_real_matrix_complex_vector_multiplication(test_report)

    call test_complex_matrix_real_vector_multiplication(test_report)

    call test_real_matrix_matrix_multiplication(test_report)

    call test_complex_matrix_matrix_multiplication(test_report)

    call test_real_matrix_complex_matrix_multiply(test_report)

    call test_complex_matrix_real_matrix_multiply(test_report)
    
    if (present(kill_on_failure)) then
      call test_report%report('general_matrix_multiplication', kill_on_failure)
    else
      call test_report%report('general_matrix_multiplication')
    end if

    call test_report%finalise()
  end subroutine general_matrix_multiplication_test_driver


  !> Test general matrix vector product for a real matrix and vector.
  subroutine test_real_matrix_vector_multiplication(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    
    real(dp), allocatable :: A(:,:), b(:), c(:), c_ref(:)

    ! setup real test arrays
    A = real_matrix_5x7
    b = real_vector_7
    c_ref = matmul(A, b)
    allocate(c(5))
    call matrix_multiply(A, b, c)
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix vector product for a real matrix and vector. &
                           Expected: A * b')
    
    b = real_vector_5
    c_ref = matmul(transpose(A), b) 
    deallocate(c); allocate(c(7))
    call matrix_multiply(A, b, c, trans_A = "T")
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix vector product for a real matrix and vector. &
                           Expected: transpose(A) * b')

  end subroutine test_real_matrix_vector_multiplication


  !> Test general matrix vector product for a complex matrix and vector.
  subroutine test_complex_matrix_vector_multiplication(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
                            
    complex(dp), allocatable :: A(:,:), b(:), c(:), c_ref(:)

    A = complex_matrix_5x7
    b = complex_vector_7
    c_ref = matmul(A, b)
    allocate(c(5))
    call matrix_multiply(A, b, c)
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix vector product for a complex matrix and vector. &
                           Expected: A * b')

    b = complex_vector_5
    c_ref = matmul(transpose(A), b)
    deallocate(c); allocate(c(7))
    call matrix_multiply(A, b, c, trans_A = "T")
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix vector product for a complex matrix and vector. &
                           Expected: transpose(A) * b')

    b = complex_vector_5
    c_ref = matmul(transpose(conjg(A)), b)
    deallocate(c); allocate(c(7))
    call matrix_multiply(A, b, c, trans_A = "C")
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix vector product for a complex matrix and vector. &
                           Expected: adjungate(A) * b')
  end subroutine test_complex_matrix_vector_multiplication


  !> Test general matrix-vector product for a real matrix and a complex vector
  subroutine test_real_matrix_complex_vector_multiplication(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    
    real(dp), allocatable :: A(:, :)
    complex(dp), allocatable :: b(:), c(:), c_ref(:)

    A = real_matrix_5x7
    b = complex_vector_7
    c_ref = matmul(A, b)
    allocate(c(5))
    call matrix_multiply(A, b, c)
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix-vector product for a real matrix A and a complex vector b. &
                           Expected: A * b')

    b = complex_vector_5
    c_ref = matmul(b, A)
    deallocate(c); allocate(c(7))
    call matrix_multiply(A, b, c, trans_A = "T")
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix-vector product for a real matrix A and complex vector b. &
                           Expected: transpose(A) * b')
  end subroutine test_real_matrix_complex_vector_multiplication


  !> Test general matrix product for a complex matrix and a real vector
  subroutine test_complex_matrix_real_vector_multiplication(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    
    real(dp), allocatable :: b(:)
    complex(dp), allocatable :: A(:, :), c(:), c_ref(:)

    A = complex_matrix_5x7
    b = real_vector_7
    allocate(c(5))
    c_ref = matmul(A, b)
    call matrix_multiply(A, b, c)
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix-vector product for a complex matrix A and a real vector b. &
                           Expected: A * b')
    
    b = real_vector_5
    deallocate(c); allocate(c(7))
    c_ref = matmul(transpose(A), b)
    call matrix_multiply(A, b, c, trans_A = "T")
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix-vector product for a complex matrix A and a real vector b. &
                           Expected: transpose(A) * b')
    
    c_ref = matmul(b, conjg(A))
    call matrix_multiply(A, b, c, trans_A = "C")
    call test_report%assert(all_close(c, c_ref), &
                           'Test general matrix-vector product for a complex matrix A and a real vector b. &
                           Expected: adjungate(A) * b')
  end subroutine test_complex_matrix_real_vector_multiplication


  subroutine test_real_matrix_matrix_multiplication(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report
                            
    real(dp), allocatable :: A(:,:), B(:,:), C(:,:), C_ref(:,:)

    A = real_matrix_5x7
    B = real_matrix_7x5
    allocate(C(5, 5))
    C_ref = matmul(A, B)
    call matrix_multiply(A, B, C)
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for real matrices. &
                           Expected: A * B')
    
    A = real_matrix_5x7
    B = real_matrix_5x7
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), B)
    call matrix_multiply(A, B, C, trans_A = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for real matrices. &
                           Expected: transpose(A) * B')

    deallocate(C); allocate(C(5, 5))
    C_ref = matmul(A, transpose(B))
    call matrix_multiply(A, B, C, trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for real matrices. &
                           Expected: A * transpose(B)')

    A = real_matrix_5x7
    B = real_matrix_7x5
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), transpose(B))
    call matrix_multiply(A, B, C, trans_A = "T", trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for real matrices. &
                           Expected: transpose(A) * transpose(B)')
  end subroutine test_real_matrix_matrix_multiplication


  !> Test general matrix-matrix product for complex matrices.
  subroutine test_complex_matrix_matrix_multiplication(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report
                            
    complex(dp), allocatable :: A(:,:), B(:,:), C(:,:), C_ref(:,:)


    ! None transposed 
    A = complex_matrix_5x7
    B = complex_matrix_7x5
    allocate(C(5, 5))
    C_ref = matmul(A, B)
    call matrix_multiply(A, B, C)
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: A * B')

    ! A transposed
    A = complex_matrix_5x7
    B = complex_matrix_5x7
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), B)

    call matrix_multiply(A, B, C, trans_A = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * B')

    C_ref = matmul(transpose(conjg(A)), B)
    call matrix_multiply(A, B, C, trans_A = "C")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: adjungate(A) * B')

    ! B transposed
    deallocate(C); allocate(C(5, 5))
    C_ref = matmul(A, transpose(B))
    call matrix_multiply(A, B, C, trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * B')

    C_ref = matmul(A, transpose(conjg(B)))
    call matrix_multiply(A, B, C, trans_B = "C")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: adjungate(A) * B')

    ! Both transposed
    A = complex_matrix_5x7
    B = complex_matrix_7x5
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), transpose(B))
    call matrix_multiply(A, B, C, trans_A = "T", trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * transpose(B)')
    

    C_ref = matmul(transpose(conjg(A)), transpose(B))
    call matrix_multiply(A, B, C, trans_A = "C", trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: ajungate(A) * transpose(B)')

    C_ref = matmul(transpose(A), transpose(conjg(B)))
    call matrix_multiply(A, B, C, trans_A = "T", trans_B = "C")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * ajungate(B)')

    C_ref = conjg(matmul(transpose(A), transpose(B)))
    call matrix_multiply(A, B, C, trans_A = "C", trans_B = "C")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: ajungate(A) * ajungate(B)')
                           
  end subroutine test_complex_matrix_matrix_multiplication


  !> Test genreal matrix-matrix product for a real matrix times a complex matrix.
  subroutine test_real_matrix_complex_matrix_multiply(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    
    real(dp), allocatable :: A(:, :)
    complex(dp), allocatable :: B(:,:), C(:,:), C_ref(:,:)

   ! None transposed 
    A = real_matrix_5x7
    B = complex_matrix_7x5
    allocate(C(5, 5))
    C_ref = matmul(A, B)
    call matrix_multiply(A, B, C)
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: A * B')

    ! A transposed
    A = real_matrix_5x7
    B = complex_matrix_5x7
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), B)
    call matrix_multiply(A, B, C, trans_A = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * B')

    ! B transposed
    deallocate(C); allocate(C(5, 5))
    C_ref = matmul(A, transpose(B))
    call matrix_multiply(A, B, C, trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * B')

    C_ref = matmul(A, conjg(transpose(B)))
    call matrix_multiply(A, B, C, trans_B = "C")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: adjungate(A) * B')

    ! Both transposed
    A = real_matrix_5x7
    B = complex_matrix_7x5
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), transpose(B))

    call matrix_multiply(A, B, C, trans_A = "T", trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * transpose(B)')

    C_ref = matmul(transpose(A), transpose(conjg(B)))
    call matrix_multiply(A, B, C, trans_A = "T", trans_B = "C")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * ajungate(B)')
  end subroutine test_real_matrix_complex_matrix_multiply


  !> Test genreal matrix-matrix product for a real matrix times a complex matrix.
  subroutine test_complex_matrix_real_matrix_multiply(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    
    real(dp), allocatable :: B(:, :)
    complex(dp), allocatable :: A(:,:), C(:,:), C_ref(:,:)

    ! None transposed 
    A = complex_matrix_5x7
    B = complex_matrix_7x5
    allocate(C(5, 5))
    C_ref = matmul(A, B)

    call matrix_multiply(A, B, C)
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: A * B')

    ! A transposed
    A = complex_matrix_5x7
    B = complex_matrix_5x7
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), B)
    call matrix_multiply(A, B, C, trans_A = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * B')

    C_ref = matmul(conjg(transpose(A)), B)
    call matrix_multiply(A, B, C, trans_A = "C") 
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: adjungate(A) * B')


    ! B transposed
    deallocate(C); allocate(C(5, 5))
    C_ref = matmul(A, transpose(B))
    call matrix_multiply(A, B, C, trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * B')

    ! Both transposed
    A = complex_matrix_5x7
    B = complex_matrix_7x5
    deallocate(C); allocate(C(7, 7))
    C_ref = matmul(transpose(A), transpose(B))
    call matrix_multiply(A, B, C, trans_A = "T", trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: transpose(A) * transpose(B)')   

    C_ref = matmul(transpose(conjg(A)), transpose(B))
    call matrix_multiply(A, B, C, trans_A = "C", trans_B = "T")
    call test_report%assert(all_close(C, C_ref), &
                           'Test general matrix-matrix product for complex matrices. &
                           Expected: ajungate(A) * transpose(B)')
  end subroutine test_complex_matrix_real_matrix_multiply

end module general_matrix_multiplication_test