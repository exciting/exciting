!> Tests for `solve_generalized_symmetric_eigenproblem`.
module generalized_symmetric_eigenproblem_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, identity_real_dp
  use mock_arrays, only: real_symmetric_matrix_5x5, real_symmetric_positive_definite_matrix_5x5
  use generalized_symmetric_eigenproblem, only: solve_generalized_symmetric_eigenproblem

  implicit none

  private
  public :: solve_generalized_symmetric_eigenproblem_test_driver

  contains

  !> Run tests for generalized eigenproblems
  subroutine solve_generalized_symmetric_eigenproblem_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 2

    call test_report%init(n_assertions, mpiglobal)

    ! Run unit tests
    call test_solve_gen_sym_eigenproblem_real_dp(test_report)

    if (present(kill_on_failure)) then
       call test_report%report('solve_generalized_symmetric_eigenproblem', kill_on_failure)
    else
       call test_report%report('solve_generalized_symmetric_eigenproblem')
    end if

    call test_report%finalise()
  end subroutine solve_generalized_symmetric_eigenproblem_test_driver

  !> Test solving a real generalized symmetric-definite eigenproblem with 5x5 matrices, where 4 eigenvalues are being
  !> searched for.
  subroutine test_solve_gen_sym_eigenproblem_real_dp(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report

    logical :: evecs_orthonormal, fullfills_eigenvalue_eq
    integer ::  N
    real(dp) :: A(5, 5), B(5, 5), eigenvalues(4), eigenvectors(5, 4)

    A = real_symmetric_matrix_5x5
    B = real_symmetric_positive_definite_matrix_5x5

    N = size(A, dim=1)

    call solve_generalized_symmetric_eigenproblem(A, B, 1e-12_dp, eigenvalues, eigenvectors)

    ! Re-set original matrices (since A and B are destroyed during the LAPACK call)
    A = real_symmetric_matrix_5x5
    B = real_symmetric_positive_definite_matrix_5x5

    ! Verify that eigenvectors are orthonormal
    evecs_orthonormal = all_close(matmul(transpose(eigenvectors), matmul(B, eigenvectors)), identity_real_dp(4), 1e-12_dp)
    call test_report%assert(evecs_orthonormal, 'Eigenvectors are not orthonormal (X^T * B * X = 1).')

    ! Verify that the eigenvalue equations A * X_i = lambda_i * B * X_i hold
    fullfills_eigenvalue_eq = all_close(matmul(A, eigenvectors), matmul(B,  spread(eigenvalues, 1, N) * eigenvectors), 1e-12_dp)
    call test_report%assert(fullfills_eigenvalue_eq, 'Eigenvalue equation is not fulfilled (A * X_i = lambda_i * B * X_i).')

  end subroutine test_solve_gen_sym_eigenproblem_real_dp

end module generalized_symmetric_eigenproblem_test
