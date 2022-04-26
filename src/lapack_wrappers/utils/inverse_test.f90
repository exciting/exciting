module inverse_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, identity_real_dp, identity_complex_dp
  use mock_arrays, only: real_full_rank_matrix_5x5, &
                         complex_full_rank_matrix_5x5
  use inverse, only: inverse_LU, invert_LU

  implicit none

  private
  public :: inverse_test_driver

contains

  !> Run tests for LU factorization
  subroutine inverse_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 8

    call test_report%init(n_assertions, mpiglobal)

    call test_inverse_LU(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('inverse', kill_on_failure)
    else
      call test_report%report('inverse')
    end if

    call test_report%finalise()
  end subroutine inverse_test_driver

  !> Test inverse
  subroutine test_inverse_LU(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: A_real(5, 5), A_inv_real(5, 5)
    complex(dp) :: A_complex(5, 5), A_inv_complex(5, 5)

    A_real = real_full_rank_matrix_5x5
    A_inv_real = inverse_LU(A_real)
    call test_report%assert(all_close(matmul(A_inv_real, A_real), identity_real_dp(5)), &
                           'Test for inverse with real input (function call) that the result is left-inverse. &
                           Expected: A^-1 * A = 1')

    call test_report%assert(all_close(matmul(A_real, A_inv_real), identity_real_dp(5)), &
                           'Test for inverse with real input (function call) that the result is right-inverse. &
                           Expected: A * A^-1 = 1')

    A_inv_real = A_real
    call invert_LU(A_inv_real)
    call test_report%assert(all_close(matmul(A_inv_real, A_real), identity_real_dp(5)), &
                           'Test for invert with real input (subroutine call) that the result is left-inverse. &
                           Expected: A^-1 * A = 1')

    call test_report%assert(all_close(matmul(A_real, A_inv_real), identity_real_dp(5)), &
                           'Test for invert with real input (subroutine call) that the result is right-inverse. &
                           Expected: A * A^-1 = 1')

    A_complex = complex_full_rank_matrix_5x5
    A_inv_complex = inverse_LU(A_complex)
    call test_report%assert(all_close(matmul(A_inv_complex, A_complex), identity_complex_dp(5)), &
                           'Test for inverse with complex input (function call) that the result is left-inverse. &
                           Expected: A^-1 * A = 1')

    call test_report%assert(all_close(matmul(A_complex, A_inv_complex), identity_complex_dp(5)), &
                           'Test for inverse with complex input (function call) that the result is right-inverse. &
                           Expected: A * A^-1 = 1')

    A_inv_complex = A_complex
    call invert_LU(A_inv_complex)
    call test_report%assert(all_close(matmul(A_inv_complex, A_complex), identity_complex_dp(5)), &
                           'Test for invert with complex input (subroutine call) that the result is left-inverse. &
                           Expected: A^-1 * A = 1')

    call test_report%assert(all_close(matmul(A_complex, A_inv_complex), identity_complex_dp(5)), &
                           'Test for invert with complex input (subroutine call) that the result is right-inverse. &
                           Expected: A * A^-1 = 1')
  end subroutine test_inverse_LU

end module inverse_test