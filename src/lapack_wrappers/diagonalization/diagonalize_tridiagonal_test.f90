!> Tests for `diagonalize_tridiagonal`.
! TODO(Bene): Issue #137: test `diagonalize_symtridiag` with `transformation_matrix_in` as optional input
module diagonalize_tridiagonal_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, identity_real_dp
  use modmpi, only: terminate_if_false
  use mock_arrays, only: real_vector_5, real_vector_7
  use diagonalize_tridiagonal, only: diagonalize_symtridiag
  
  implicit none

  private
  public :: diagonalize_tridiagonal_test_driver
  
contains

  !> Run tests for tridiagonal T diagonalization
  subroutine diagonalize_tridiagonal_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 3

    call test_report%init(n_assertions, mpiglobal)

    ! Run unit tests

    call test_diagonalize_real_symmetric_tridiagonal(test_report)
      
    if (present(kill_on_failure)) then
      call test_report%report('diagonalize_tridiagonal', kill_on_failure)
    else
      call test_report%report('diagonalize_tridiagonal')
    end if

    call test_report%finalise()
  end subroutine diagonalize_tridiagonal_test_driver

  !> Test diagonalization for a real symmetric tridiagonal T.
  subroutine test_diagonalize_real_symmetric_tridiagonal(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report 

    logical :: evecs_orthornormal, fullfills_eigenvalue_eq
    integer :: i, N
    real(dp), allocatable :: lambda(:), lambda_only(:), X(:, :), T(:, :)

    real(dp), parameter :: diagonal(6) = [-26.54377445_dp, -17.22046184_dp,   4.83382119_dp,  -9.47731147_dp, &
            0.03664525_dp,   3.04736258_dp]
    real(dp), parameter :: subdiagonal(5) = [ 31.18072641_dp,  -4.72689264_dp,  48.67919055_dp, -41.66440826_dp, &
            -40.21374565_dp]


    N = size(diagonal)

    call diagonalize_symtridiag(diagonal, subdiagonal, lambda_only)
    call diagonalize_symtridiag(diagonal, subdiagonal, lambda, X)

    call test_report%assert(all_close(lambda_only, lambda), 'diagonalize_symtridiag does not return the same eigen &
                           values whith and without calculating the eigen vectors.')

    evecs_orthornormal = all_close(matmul(transpose(X), X), identity_real_dp(size(diagonal)))
    call test_report%assert(evecs_orthornormal, 'Eigen vectors returned by diagonalize_symtridiag are not orthornormal &
                           (X^T * X = 1).')
    
    allocate(T(N, N))
    T = 0.0_dp
    T(1, 1) = diagonal(1)
    do i=2, N
      T(i, i) = diagonal(i)
      T(i, i-1) = subdiagonal(i-1)
      T(i-1, i) = subdiagonal(i-1)
    end do 

    fullfills_eigenvalue_eq = all_close(matmul(T, X), spread(lambda, 1, N) * X)
    call test_report%assert(fullfills_eigenvalue_eq, 'Eigen value equation is not fulfilled for results of &
            diagonalize_symtridiag (T * x_i = lambda_i * x_i).')
  end subroutine test_diagonalize_real_symmetric_tridiagonal

end module diagonalize_tridiagonal_test
