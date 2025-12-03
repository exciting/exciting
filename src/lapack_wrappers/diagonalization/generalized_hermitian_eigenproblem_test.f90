!> Tests for `solve_generalized_symmetric_eigenproblem`.
module generalized_hermitian_eigenproblem_test
  use generalized_hermitian_eigenproblem, only: solve_generalized_hermitian_eigenproblem
  use math_utils, only: all_close, identity_real_dp, identity_complex_dp
  use mock_arrays, only: real_symmetric_matrix_5x5, real_positive_definite_matrix_5x5, &
                         complex_hermitian_matrix_5x5, complex_positive_definite_matrix_5x5
  use modmpi, only: mpiinfo
  use precision, only: dp, i32
  use unit_test_framework, only: unit_test_type
  use to_char_conversion, only: to_char

  implicit none

  private
  public :: solve_generalized_hermitian_eigenproblem_test_driver

  contains

  !> Run tests for generalized hermitian eigenproblems
  subroutine solve_generalized_hermitian_eigenproblem_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test object
    type(unit_test_type) :: test_report

    call test_report%init( mpiglobal)

    ! Run unit tests
    call test_solve_gen_hermitian_eigenproblem_complex_dp(test_report)
    call test_solve_gen_sym_eigenproblem_real_dp(test_report)
    
    if (present(kill_on_failure)) then
       call test_report%report('solve_generalized_symmetric_eigenproblem', kill_on_failure)
    else
       call test_report%report('solve_generalized_symmetric_eigenproblem')
    end if

    call test_report%finalise()
  end subroutine solve_generalized_hermitian_eigenproblem_test_driver

  !> Test solving a real generalized symmetric-definite eigenproblem with 5x5 matrices, where 4 eigenvalues are being
  !> searched for.
  subroutine test_solve_gen_sym_eigenproblem_real_dp(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report

    real(dp), parameter :: tol = 1e-12_dp
    logical :: evecs_orthonormal, fullfills_eigenvalue_eq
    integer ::  N
    real(dp) :: A(5, 5), B(5, 5), eigenvalues(4), eigenvectors(5, 4)

    A = real_symmetric_matrix_5x5
    B = real_positive_definite_matrix_5x5

    N = size(A, dim=1)

    call solve_generalized_hermitian_eigenproblem(A, B, tol, eigenvalues, eigenvectors)

    ! Re-set original matrices (since A and B are destroyed during the LAPACK call)
    A = real_symmetric_matrix_5x5
    B = real_positive_definite_matrix_5x5

    ! Verify that eigenvectors are orthonormal
    evecs_orthonormal = all_close(matmul(transpose(eigenvectors), matmul(B, eigenvectors)), identity_real_dp(4), tol)
    call test_report%assert(evecs_orthonormal, 'Eigenvectors are not orthonormal (X^T * B * X = 1).')

    ! Verify that the eigenvalue equations A * X_i = lambda_i * B * X_i hold
    fullfills_eigenvalue_eq = all_close(matmul(A, eigenvectors), matmul(B,  spread(eigenvalues, 1, N) * eigenvectors), tol)
    call test_report%assert(fullfills_eigenvalue_eq, '(symmetric_eigenproblem) Eigenvalue equation is not fulfilled (A * X_i = lambda_i * B * X_i).')

  end subroutine test_solve_gen_sym_eigenproblem_real_dp


  !> Test solving a complex generalized hermitian eigenproblem with 5x5 matrices, 
  !> where 4 eigenvalues are being searched for.
  subroutine test_solve_gen_hermitian_eigenproblem_complex_dp(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report

    real(dp), parameter :: tol_solver = 1e-14_dp
    real(dp), parameter :: tol_test = 1e-11_dp
    logical :: evecs_orthonormal, fullfills_eigenvalue_eq
    integer(i32) ::  N
    complex(dp) :: A(5, 5), B(5, 5), eigenvectors(5, 4)
    real(dp) :: eigenvalues(4)
    real(dp), allocatable :: eigenvalues_ref(:)

    A = complex_hermitian_matrix_5x5
    B = complex_positive_definite_matrix_5x5
    N = size(A, dim=1)

    call solve_generalized_hermitian_eigenproblem(A, B, tol_solver, eigenvalues, eigenvectors)

    ! Re-set original matrices (since A and B are destroyed during the LAPACK call)
    A = complex_hermitian_matrix_5x5
    B = complex_positive_definite_matrix_5x5

    ! Verify that eigenvectors are orthonormal
    evecs_orthonormal = all_close(matmul(conjg(transpose(eigenvectors)), matmul(B, eigenvectors)), identity_complex_dp(4), tol_test)
    call test_report%assert(evecs_orthonormal, 'Eigenvectors are not orthonormal (X^H * B * X = 1).')

    ! Verify that the eigenvalue equations A * X_i = lambda_i * B * X_i hold
    fullfills_eigenvalue_eq = all_close(matmul(A, eigenvectors), matmul(B,  spread(eigenvalues, 1, N) * eigenvectors), tol_test)
    call test_report%assert(fullfills_eigenvalue_eq, '(hermitian_eigenproblem) Eigenvalue equation is not fulfilled (A * X_i = lambda_i * B * X_i). Diff = ' // &
      to_char( maxval( abs(matmul(A, eigenvectors) - matmul(B,  spread(eigenvalues, 1, N) * eigenvectors)) ) ) )

    ! Obtain only the eigenvalues
    eigenvalues_ref = eigenvalues
    call solve_generalized_hermitian_eigenproblem(A, B, tol_solver, eigenvalues)
    call test_report%assert( all_close(eigenvalues, eigenvalues_ref, tol=tol_test), '(hermitian_eigenproblem) Eigenvalues are not equal to the reference')

  end subroutine

end module generalized_hermitian_eigenproblem_test
