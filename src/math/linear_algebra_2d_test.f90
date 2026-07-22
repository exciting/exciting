module linear_algebra_2d_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use math_utils, only: all_close, all_zero
  use linear_algebra_2d, only: solve_2d_cramer

  implicit none

  private
  public :: linear_algebra_2d_test_driver
    
  
  contains

  !> Run tests for math tools
  subroutine linear_algebra_2d_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure  
    !> test object
    type(unit_test_type) :: test_report

    ! Initialize test object
    call test_report%init(mpiglobal)

    ! Run and assert tests
    
    call test_solve_2d_cramer(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('linear_algebra_2d', kill_on_failure)
    else
      call test_report%report('linear_algebra_2d')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine linear_algebra_2d_test_driver


  !> Test solve_2x2_cramer
  subroutine test_solve_2d_cramer(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    !> test matrices and vectors
    real(dp) :: A(2,2), b(2), x(2), x_expected(2)
    real(dp) :: A_singular(2,2), b_singular(2)
    integer :: i, info
    real(dp), parameter :: rtol = 1e-12_dp
    real(dp) :: atol

    ! -------------------------------------------------------------
    ! Case 1: simple system with exact integer solution
    ! [ 2  1 ] [x1] = [5]
    ! [ 1  3 ] [x2]   [6]
    ! Solution: x1 = 1.8, x2 = 1.4
    ! -------------------------------------------------------------
    A = reshape([2.0_dp, 1.0_dp, &
                 1.0_dp, 3.0_dp], [2,2])
    b = [5.0_dp, 6.0_dp]
    x_expected = [1.8_dp, 1.4_dp]

    call solve_2d_cramer(A, b, x, info)

    call test_report%assert(all_close(x, x_expected, 1.0e-12_dp), &
        'Test solve_2d_cramer for small real system. Expected x = [1.8, 1.4].')

    ! -------------------------------------------------------------
    ! Case 2: identity matrix — solution should be equal to b
    ! -------------------------------------------------------------
    A = reshape([1.0_dp, 0.0_dp, &
                 0.0_dp, 1.0_dp], [2,2])
    b = [3.14_dp, -2.71_dp]

    call solve_2d_cramer(A, b, x, info)

    call test_report%assert(all_close(x, b), &
        'Test solve_2d_cramer with identity matrix. Expected x = b.')

    ! -------------------------------------------------------------
    ! Case 3: scaling test — matrix with large and small values
    ! -------------------------------------------------------------
    A = reshape([1.0e-3_dp, 2.0e-3_dp, &
                 4.0e+3_dp, 5.0e+3_dp], [2,2])
    b = [7.0e-3_dp, 8.0e+3_dp]

    atol = rtol * maxval(abs(b))

    call solve_2d_cramer(A, b, x, info)
    ! Verify by re-multiplying
    call test_report%assert(all_close(matmul(A, x), b, atol), &
        'Test solve_2d_cramer numerical consistency (A*x ≈ b).')

    ! -------------------------------------------------------------
    ! Case 4: singular matrix — determinant = 0 should raise an error
    ! -------------------------------------------------------------
    A = reshape([1.0_dp, 2.0_dp, &
                          2.0_dp, 4.0_dp], [2,2])
    b = [1.0_dp, 2.0_dp]

    call solve_2d_cramer(A, b, x, info)
    ! Verify that an error is raised
    call test_report%assert(info==1, &
        'Test solve_2d_cramer determinant=0.')

    ! -------------------------------------------------------------
    ! Case 5: both A and b are very small (order 1e-50)
    ! -------------------------------------------------------------
    A = reshape([1.0e-50_dp, 2.0e-50_dp, &
                 3.0e-50_dp, 4.0e-50_dp], [2,2])
    b = [5.0e-50_dp, 11.0e-50_dp]

    atol = rtol * maxval(abs(b))

    call solve_2d_cramer(A, b, x, info)
    call test_report%assert(all_close(matmul(A, x), b, atol), &
        'Test solve_2d_cramer with very small A and b (order 1e-50).')

    ! -------------------------------------------------------------
    ! Case 6: both A and b are very large (order 1e+50)
    ! -------------------------------------------------------------
    A = reshape([1.0e+50_dp, 2.0e+50_dp, &
                 3.0e+50_dp, 4.0e+50_dp], [2,2])
    b = [5.0e+50_dp, 11.0e+50_dp]

    atol = rtol * maxval(abs(b))

    call solve_2d_cramer(A, b, x, info)
    call test_report%assert(all_close(matmul(A, x), b, atol), &
        'Test solve_2d_cramer with very large A and b (order 1e+50).')

    ! -------------------------------------------------------------
    ! Case 7: A is very large, b is very small
    ! -------------------------------------------------------------
    A = reshape([1.0e+50_dp, 2.0e+50_dp, &
                 3.0e+50_dp, 4.0e+50_dp], [2,2])
    b = [5.0e-50_dp, 11.0e-50_dp]

    atol = rtol * maxval(abs(b))

    call solve_2d_cramer(A, b, x, info)
    call test_report%assert(all_close(matmul(A, x), b, atol),&
        'Test solve_2d_cramer with large A and small b.')

    ! -------------------------------------------------------------
    ! Case 8: A is very small, b is very large
    ! -------------------------------------------------------------
    A = reshape([1.0e-50_dp, 2.0e-50_dp, &
                 3.0e-50_dp, 4.0e-50_dp], [2,2])
    b = [5.0e+50_dp, 11.0e+50_dp]

    call solve_2d_cramer(A, b, x, info)

    atol = rtol * maxval(abs(b))

    call test_report%assert(all_close(matmul(A, x), b, atol),&
        'Test solve_2d_cramer with small A and large b.')


  end subroutine test_solve_2d_cramer


end module linear_algebra_2d_test