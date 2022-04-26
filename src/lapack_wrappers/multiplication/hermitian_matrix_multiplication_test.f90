module hermitian_matrix_multiplication_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close
  use mock_arrays, only: real_vector_5, &
                         complex_vector_5, &
                         real_matrix_5x7, &
                         real_matrix_7x5, &
                         real_symmetric_matrix_5x5, &
                         complex_hermitian_matrix_5x5, &
                         complex_matrix_5x7, &
                         complex_matrix_7x5
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply

  implicit none

  private
  public :: hermitian_matrix_multiplication_test_driver

contains

  !> Run tests for hermitian matrix multiplication
  subroutine hermitian_matrix_multiplication_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 24

    call test_report%init(n_assertions, mpiglobal)

    call test_real_symmetric_matrix_vector_multiplication(test_report)

    call test_complex_hermitian_matrix_vector_multiplication(test_report)

    call test_real_symmetric_matrix_complex_vector_multiplication(test_report)

    call test_complex_hermitian_matrix_real_vector_multiplication(test_report)

    call test_real_symmetric_matrix_matrix_multiplication(test_report)

    call test_complex_hermitian_matrix_matrix_multiplication(test_report)

    call test_real_symmetric_matrix_complex_matrix_multiplication(test_report)

    call test_complex_hermitian_matrix_real_matrix_multiplication(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('hermitian_matrix_multiplication', kill_on_failure)
    else
      call test_report%report('hermitian_matrix_multiplication')
    end if

    call test_report%finalise()
  end subroutine hermitian_matrix_multiplication_test_driver


  !> Test hermitian matrix-vector product for a real symmetric matrix and a real vector. 
  subroutine test_real_symmetric_matrix_vector_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    real(dp), allocatable :: A(:, :), b(:), c(:), c_ref(:)

    A = real_symmetric_matrix_5x5
    b = real_vector_5
    c_ref = matmul(A, b)
    allocate(c(5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, b, c)
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix-vector product for a real matrix and vector &
                          where the upper triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'L')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix-vector product for a real matrix and vector &
                          where the upper triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')
  end subroutine test_real_symmetric_matrix_vector_multiplication


  !> Test hermitian matrix-vector product for a complex hermitian matrix and a complex vector.
  subroutine test_complex_hermitian_matrix_vector_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    complex(dp), allocatable :: A(:, :), b(:), c(:), c_ref(:)

    integer :: i, j

    A = complex_hermitian_matrix_5x5
    b = complex_vector_5
    c_ref = matmul(A, b)
    allocate(c(5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'U')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix-vector product for a complex hermitian matrix and a complex vector &
                          where the upper triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'L')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix-vector product for a complex hermitian matrix and a complex vector &
                          where the lower triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')
  end subroutine test_complex_hermitian_matrix_vector_multiplication


  !> Test hermitian matrix-vector product for a real symmetric matrix and a complex vector.
  subroutine test_real_symmetric_matrix_complex_vector_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    real(dp), allocatable :: A(:, :)
    complex(dp), allocatable :: b(:), c(:), c_ref(:)

    A = real_symmetric_matrix_5x5
    b = complex_vector_5
    c_ref = matmul(A, b)
    allocate(c(5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'u')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix-vector product for a real symmetric matrix and a complex vector &
                          where the upper triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'L')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix-vector product for a real symmetric matrix and a complex vector &
                          where the lower triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')
  end subroutine test_real_symmetric_matrix_complex_vector_multiplication


  !> Test hermitian matrix-vector product for a complex hermitian matrix and a real vector.
  subroutine test_complex_hermitian_matrix_real_vector_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    real(dp), allocatable :: b(:)
    complex(dp), allocatable :: A(:, :), c(:), c_ref(:)

    A = complex_hermitian_matrix_5x5
    b = real_vector_5
    c_ref = matmul(A, b)
    allocate(c(5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'U')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix vector product for a complex hermitian matrix and a real vector &
                          where the upper triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, b, c, uplo = 'L')
    call test_report%assert(all_close(c, c_ref), &
                          'Test hermitian matrix vector product for a complex hermitian matrix and a real vector &
                          where the lower triangular part of the matrix is refered. &
                          Expected: Same result as general matrix-vector multiplication')
  end subroutine test_complex_hermitian_matrix_real_vector_multiplication
  

  !> Test hermitian matrix-matrix product for a real symmetric matrix and a real vector.
  subroutine test_real_symmetric_matrix_matrix_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    real(dp), allocatable :: A(:, :), B(:, :), C(:, :), C_ref(: ,:)

    ! A is symmetric
    A = real_symmetric_matrix_5x5
    B = real_matrix_5x7
    C_ref = matmul(A, B)
    allocate(C(5, 7))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'U', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a real symmetric matrix and a real vector &
                           where the left matrix is symmetric &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'L', side='l')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a real symmetric matrix and a real vector &
                           where the left matrix is symmetric &
                           and the lower triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! B is symmetric
    A = real_matrix_7x5
    B = real_symmetric_matrix_5x5
    C_ref = matmul(A, B)
    deallocate(C); allocate(C(7, 5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'u', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a real symmetric matrix and a real vector &
                           where the right matrix is symmetric &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'l', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a real symmetric matrix and a real vector &
                           where the right matrix is symmetric &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')
  end subroutine test_real_symmetric_matrix_matrix_multiplication


  !> Test hermitian matrix-matrix product for a complex hermitian matrix and a complex matrix.
  subroutine test_complex_hermitian_matrix_matrix_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    complex(dp), allocatable :: A(:, :), B(:, :), C(:, :), C_ref(: ,:)

    ! A is hermitian
    A = complex_hermitian_matrix_5x5
    B = complex_matrix_5x7
    C_ref = matmul(A, B)
    allocate(C(5, 7))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'U', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'l', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! B is hermitian
    A = complex_matrix_7x5
    B = complex_hermitian_matrix_5x5
    C_ref = matmul(A, B)
    deallocate(C); allocate(C(7, 5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'u', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'L', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')
  end subroutine test_complex_hermitian_matrix_matrix_multiplication


  !> Test hermitian matrix product for a real symmetric matrix and a complex matrix. 
  subroutine test_real_symmetric_matrix_complex_matrix_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    real(dp), allocatable :: A(:, :)
    complex(dp), allocatable :: B(:, :), C(:, :), C_ref(: ,:)

    ! A is symmetric
    A = real_symmetric_matrix_5x5
    B = complex_matrix_5x7
    C_ref = matmul(A, B)
    allocate(C(5, 7))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'U', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'l', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! B is symmetric
    A = complex_matrix_7x5
    B = real_symmetric_matrix_5x5
    C_ref = matmul(A, B)
    deallocate(C); allocate(C(7, 5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'u', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'L', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')
  end subroutine test_real_symmetric_matrix_complex_matrix_multiplication


  !> Test hermitian matrix product for a complex hermitian matrix and a real matrix.
  subroutine test_complex_hermitian_matrix_real_matrix_multiplication(test_report)
    !> Test report
    type(unit_test_type) :: test_report
  
    real(dp), allocatable :: B(:, :)
    complex(dp), allocatable :: A(:, :), C(:, :), C_ref(: ,:)

    ! A is hermitian
    A = complex_hermitian_matrix_5x5
    B = real_matrix_5x7
    C_ref = matmul(A, B)
    allocate(C(5, 7))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'U', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'l', side='L')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! B is hermitian
    A = real_matrix_7x5
    B = complex_hermitian_matrix_5x5
    C_ref = matmul(A, B)
    deallocate(C); allocate(C(7, 5))

    ! Use upper triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'u', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')

    ! Use lower triangular part
    call hermitian_matrix_multiply(A, B, C, uplo = 'L', side='R')
    call test_report%assert(all_close(C, C_ref), &
                           'Test hermitian matrix-matrix product for a complex hermitian matrix and a hermitian vector &
                           where the left matrix is hermitian &
                           and the upper triangular part of it is refered. &
                           Expected: Same result as general matrix-matrix product.')
  end subroutine test_complex_hermitian_matrix_real_matrix_multiplication

end module hermitian_matrix_multiplication_test