!> Unit tests for matrix decomposition
module lu_factorization_test
  use precision, only: dp
  use constants, only: zzero
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, is_square, identity_real_dp, identity_complex_dp
  use mock_arrays, only: real_full_rank_matrix_5x5, &
                         real_matrix_5x7, &
                         real_matrix_7x5, &
                         complex_full_rank_matrix_5x5, &
                         complex_matrix_5x7, &
                         complex_matrix_7x5
  use lu_factorization, only: lu_full_pivot, lu_row_pivot, pivot_to_permutation

  implicit none

  private
  public :: lu_factorization_test_driver
  
contains

  !> Run tests for LU factorization
  subroutine lu_factorization_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 9

    call test_report%init(n_assertions, mpiglobal)

    call test_lu_full_pivot_real(test_report)
    
    call test_lu_full_pivot_complex(test_report)

    call test_lu_row_pivot_real(test_report)

    call test_lu_row_pivot_complex(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('lu_factorization', kill_on_failure)
    else
      call test_report%report('lu_factorization')
    end if

    call test_report%finalise()
  end subroutine lu_factorization_test_driver


   !> Test LU factorization with full pivoting for real matrices
  subroutine test_lu_full_pivot_real(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    real(dp) :: A(5, 5), L(5, 5), U(5, 5)
    integer :: P(5), Q(5)

    A = real_full_rank_matrix_5x5
    call lu_full_pivot(A, L, U, P, Q, calculate_permutation=.true.)
    
    call test_report%assert(all_close(A(P, Q), matmul(L, U)), &
                           'Test for lu_full_pivot that the definition is fullfilled (real input). &
                           Expected: P * A * Q = L * U ')
  end subroutine test_lu_full_pivot_real


  !> Test LU factorization with full pivoting for complex matrices
  subroutine test_lu_full_pivot_complex(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    complex(dp) :: A(5, 5), L(5, 5), U(5, 5)
    integer :: P(5), Q(5)

    A = complex_full_rank_matrix_5x5
    call lu_full_pivot(A, L, U, P, Q, calculate_permutation=.true.)

    call test_report%assert(all_close(A(P, Q), matmul(L, U)), &
                           'Test for lu_full_pivot that the definition is fullfilled (complex input). &
                           Expected: P * A * Q = L * U ')

  end subroutine test_lu_full_pivot_complex


  !> Test LU factorization with row pivoting for real matrices
  subroutine test_lu_row_pivot_real(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    real(dp), allocatable :: A(:, :), L(:, :), U(:, :)
    integer, allocatable :: P(:)

    A = real_full_rank_matrix_5x5
    allocate(L(5, 5)); allocate(U(5, 5)); allocate(P(5))
    call lu_row_pivot(A, L, U, P, calculate_permutation=.true.)
    
    call test_report%assert(all_close(A(P, :), matmul(L, U)), &
                           'Test for lu_row_pivot that the definition is fullfilled for a real quadratic matrix. &
                           Expected: P * A = L * U ')

    deallocate(A, U)
    A = real_matrix_5x7
    allocate(U(5, 7))
    call lu_row_pivot(A, L, U, P, calculate_permutation=.true.)

    call test_report%assert(all_close(A(P, :), matmul(L, U)), &
                           'Test for lu_row_pivot that the definition is fullfilled for real matrix with m < n. &
                           Expected: P * A = L * U ')

    deallocate(A, U, L, P)
    A = real_matrix_7x5
    allocate(L(7, 5)); allocate(U(5, 5)); allocate(P(7))
    call lu_row_pivot(A, L, U, P, calculate_permutation=.true.)

    call test_report%assert(all_close(A(P, :), matmul(L, U)), &
                           'Test for lu_row_pivot that the definition is fullfilled for real matrix with m > n. &
                           Expected: P * A = L * U ')
  end subroutine test_lu_row_pivot_real

  !> Test LU factorization with row pivoting for complex matrices
  subroutine test_lu_row_pivot_complex(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    complex(dp), allocatable :: A(:, :), L(:, :), U(:, :)
    integer, allocatable :: P(:)

    A = complex_full_rank_matrix_5x5
    allocate(L(5, 5)); allocate(U(5, 5)); allocate(P(5))
    call lu_row_pivot(A, L, U, P, calculate_permutation=.true.)

    call test_report%assert(all_close(A(P, :), matmul(L, U)), &
                           'Test for lu_row_pivot that the definition is fullfilled for a complex quadratic matrix. &
                           Expected: P * A = L * U ')

    deallocate(A, U)
    A = complex_matrix_5x7
    allocate(U(5, 7))
    call lu_row_pivot(A, L, U, P, calculate_permutation=.true.)

    call test_report%assert(all_close(A(P, :), matmul(L, U)), &
                           'Test for lu_row_pivot that the definition is fullfilled for complex matrix with m < n. &
                           Expected: P * A = L * U ')

    deallocate(A, U, L, P)
    A = complex_matrix_7x5
    allocate(L(7, 5)); allocate(U(5, 5));  allocate(P(7))
    call lu_row_pivot(A, L, U, P, calculate_permutation=.true.)
    
    call test_report%assert(all_close(A(P, :), matmul(L, U)), &
                           'Test for lu_row_pivot that the definition is fullfilled for complex matrix with m > n. &
                           Expected: P * A = L * U ')
  end subroutine test_lu_row_pivot_complex

  subroutine test_pivot_to_permutation(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    integer :: pivot_map(4), permutation_map(4)

    pivot_map = [2, 2, 4, 1]
    call pivot_to_permutation(pivot_map, permutation_map)

    call test_report%assert(all(permutation_map == [4, 2, 3, 1]), &
                            'Test pivot_to_permutation for the example from the documentation. &
                            Expected: [2, 2, 4, 1] --> [4, 2, 3, 1]')

  end subroutine test_pivot_to_permutation
  
end module lu_factorization_test