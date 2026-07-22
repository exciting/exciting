module iterative_solver_test
  use precision, only: dp, i32
  use constants, only: zone, zzero
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close
  use modmpi, only: mpiinfo
#include "asserts.fpp"
  use iterative_solver, only: lanczos
  use xlapack, only: diagonalize_symtridiag

  private
  public :: iterative_solver_test_driver

  contains

  !> Run tests for iterative solver module
  subroutine iterative_solver_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> Test report object
    type(unit_test_type) :: test_report

    character(*), parameter :: test_name = 'iterative_solver'

    ! Initialize test object
    call test_report%init(mpiglobal)

    ! Run and assert tests
    call test_lanczos(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('iterative_solver', kill_on_failure)
    else
      call test_report%report('iterative_solver')
    end if
  end subroutine iterative_solver_test_driver

  !> Test the lanczos solver
  subroutine test_lanczos(test_report)
    !> Test report
    type(unit_test_type) :: test_report

    real(dp), allocatable :: alpha(:), beta(:), evals(:), evecs(:, :)
    complex(dp), allocatable :: Q_k(:, :), q_1(:)

    ! Zero tolerance
    real(dp), parameter :: zero_tol = 5e-8

    ! Expected results coming from a test run with a verified version.
    real(dp), parameter :: alpha_ref(4) = [3.0000000000000000, 3.0000000000000000, 3.0000000000000009, &
            2.9999999999999996]
    real(dp), parameter :: beta_ref(4) = [1.4142135623730951, 1.1832159566199232, 1.0141851056742199, &
            0.75592894601845451 ]
    complex(dp), parameter :: Q_k_ref(5, 4) = reshape( &
           [0.44721359549995793,  0.44721359549995793,  0.44721359549995793,  0.44721359549995793, &
            0.44721359549995793, -0.63245553203367588, -0.31622776601683794,  0.00000000000000000, &
            0.31622776601683789,  0.63245553203367588,  0.53452248382484868, -0.26726124191242445, &
           -0.53452248382484879, -0.26726124191242445,  0.53452248382484890, -0.31622776601683822, &
            0.63245553203367622,  0.00000000000000000, -0.63245553203367555,  0.31622776601683761], [5, 4]) * zone

    ! Test break down of lanczos in the first iteration
    q_1 = [zone, zzero, zzero]
    call lanczos(3, mock_zero, q_1, alpha, beta)
    call test_report%assert(.not. allocated(alpha), &
            'alpha is allocated for Lanczos breaking down in first iteration.')
    call test_report%assert(.not. allocated(beta), &
            'beta is allocated for Lanczos breaking down in first iteration.')

    call lanczos(3, mock_zero, q_1, alpha, beta, Q_k)
    call test_report%assert(.not. allocated(alpha), &
            'alpha is allocated for Lanczos breaking down in first iteration.')
    call test_report%assert(.not. allocated(beta), &
            'beta is allocated for Lanczos breaking down in first iteration.')
    call test_report%assert(.not. allocated(Q_k), &
            'Q_k is allocated for Lanczos breaking down in first iteration.')

    ! Test successful Lanczos iteration
    q_1 = zone * [1, 1, 1, 1, 1]
    call lanczos(5, mock_diagonal, q_1, alpha, beta)
    call test_report%assert(all_close(alpha, alpha_ref, zero_tol), &
            'alpha is not the same as the reference value.')
    call test_report%assert(all_close(beta, beta_ref, zero_tol), &
            'beta is not the same as the reference value.')

    call lanczos(5, mock_diagonal, q_1, alpha, beta, Q_k)
    call test_report%assert(all_close(alpha, alpha_ref, zero_tol), &
            'alpha is not the same as the reference value.')
    call test_report%assert(all_close(beta, beta_ref, zero_tol), &
            'beta is not the same as the reference value.')
    call test_report%assert(all_close(Q_k, Q_k_ref, zero_tol), &
            'Q_k is not the same as the reference value.')

    contains

    !> Mock matrix vector product that maps every vector to zero to test gringe case of
    !> lanczos breaking down after a single iteration.
    subroutine mock_zero(v_in, v_out)
      complex(dp), intent(in) :: v_in(:)
      complex(dp), intent(out) :: v_out(:)
      v_out = zzero
    end subroutine mock_zero

    !> Mock diagonal multiplication
    subroutine mock_diagonal(v_in, v_out)
      complex(dp), intent(in) :: v_in(:)
      complex(dp), intent(out) :: v_out(:)

      complex(dp), parameter :: diagonal(5) = zone * [1, 2, 3, 4, 5]

      CALL_ASSERT(size(v_in) == 5)
      CALL_ASSERT(size(v_out) == 5)

      v_out = v_in * diagonal
    end subroutine mock_diagonal

  end subroutine test_lanczos
end module iterative_solver_test
