module qr_factorization_test
  use precision, only: dp
  use constants, only:
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, all_zero, &
                        identity_real_dp, identity_complex_dp
  use mock_arrays, only: real_matrix_5x5, &
                         real_matrix_5x7, &
                         real_matrix_7x5, &
                         complex_matrix_5x5, &
                         complex_matrix_5x7, &
                         complex_matrix_7x5
  use qr_factorization, only: qr_column_pivot

  implicit none

  private
  public :: qr_factorization_test_driver

  contains

  !> Run tests for QR factorization
  subroutine qr_factorization_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 12

    call test_report%init(n_assertions, mpiglobal)

    call test_qr_column_pivot_real(test_report)
   
    call test_qr_column_pivot_complex(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('qr_factorization', kill_on_failure)
    else
      call test_report%report('qr_factorization')
    end if

    call test_report%finalise()
  end subroutine qr_factorization_test_driver
 

  !> Test QR factorization with column pivot for real input
  subroutine test_qr_column_pivot_real(test_report)
    !> Test report object
    type(unit_test_type) :: test_report
    real(dp), allocatable :: A(:, :), Q(:, :), R(:, :), L(:, :)
    integer, allocatable :: P(:)

    ! Test QR factorization for a square matrix
    A = real_matrix_5x5
    allocate(P(5)); allocate(Q(5, 5)); allocate(R(5, 5))
    call qr_column_pivot(A, P, Q, R)

    call test_report%assert(all_close(matmul(transpose(Q), Q), identity_real_dp(5)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for real input with N_rows = N_cols. &
                           Expected: Q**T * Q == I')

    call test_report%assert(all_close(matmul(Q, transpose(Q)), identity_real_dp(5)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for real input with N_rows = N_cols. &
                           Expected: Q * Q**T == I')

    call test_report%assert(all_close(A(:, P), matmul(Q, R)), &
                           'Test if qr_column_pivot for real input with N_rows = N_cols &
                           calculates Q, R and P correctly. &
                           Expected: A (:, P) == Q * R.')
    deallocate(Q, R)


    ! Test QR factorization for a matrix with N_rows > N_cols
    A = real_matrix_7x5
    allocate(Q(7, 7)); allocate(R(7, 5))
    call qr_column_pivot(A, P, Q, R)

    call test_report%assert(all_close(matmul(transpose(Q), Q), identity_real_dp(7)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for real input with N_rows > N_cols. &
                           Expected: Q**T * Q == I')

    call test_report%assert(all_close(matmul(Q, transpose(Q)), identity_real_dp(7)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for real input with N_rows > N_cols. &
                           Expected: Q * Q**T == I')

    call test_report%assert(all_close(A(:, P), matmul(Q, R)), &
                           'Test if qr_column_pivot for real input with N_rows > N_cols &
                           calculates Q, R and P correctly. &
                           Expected: A (:, P) == Q * R.')

  end subroutine test_qr_column_pivot_real


  !> Test QR factorization with column pivot for complex input
  subroutine test_qr_column_pivot_complex(test_report)
    !> Test report object
    type(unit_test_type) :: test_report
    complex(dp), allocatable :: A(:, :), Q(:, :), R(:, :), L(:, :)
    integer, allocatable :: P(:)

    ! Test QR factorization for a square matrix
    A = complex_matrix_5x5
    allocate(P(5)); allocate(Q(5, 5)); allocate(R(5, 5))
    call qr_column_pivot(A, P, Q, R)

    call test_report%assert(all_close(matmul(transpose(conjg(Q)), Q), identity_complex_dp(5)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for complex input with N_rows = N_cols. &
                           Expected: Q**T * Q == I')

    call test_report%assert(all_close(matmul(Q, transpose(conjg(Q))), identity_complex_dp(5)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for complex input with N_rows = N_cols. &
                           Expected: Q * Q**T == I')

    call test_report%assert(all_close(A(:, P), matmul(Q, R)), &
                           'Test if qr_column_pivot for complex input with N_rows = N_cols &
                           calculates Q, R and P correctly. &
                           Expected: A (:, P) == Q * R.')
    deallocate(Q, R)

    ! Test QR factorization for a matrix with N_rows > N_cols
    A = complex_matrix_7x5
    allocate(Q(7, 7)); allocate(R(7, 5))
    call qr_column_pivot(A, P, Q, R)

    call test_report%assert(all_close(matmul(transpose(conjg(Q)), Q), identity_complex_dp(7)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for complex input with N_rows > N_cols. &
                           Expected: Q**T * Q == I')

    call test_report%assert(all_close(matmul(Q, transpose(conjg(Q))), identity_complex_dp(7)), &
                           'Test if Q returned from qr_column_pivot is unitary &
                           for complex input with N_rows > N_cols. &
                           Expected: Q * Q**T == I')

    call test_report%assert(all_close(A(:, P), matmul(Q, R)), &
                           'Test if qr_column_pivot for complex input with N_rows > N_cols &
                           calculates Q, R and P correctly. &
                           Expected: A (:, P) == Q * R.')

  end subroutine test_qr_column_pivot_complex

end module qr_factorization_test