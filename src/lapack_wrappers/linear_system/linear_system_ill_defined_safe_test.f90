!> Unit tests for linear_system_ill_defined_safe module
module linear_system_ill_defined_safe_test
  use constants, only: zone, zzero
  use exciting_mpi, only: mpiinfo
  use linear_system_ill_defined_safe, only: ill_defined_safe_solve
  use math_utils, only: all_close
  use mock_arrays, only: complex_positive_definite_matrix_5x5, complex_matrix_5x5, complex_matrix_5x7
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
  use unit_test_framework, only: unit_test_type

  implicit none

  private 

  public :: linear_system_ill_defined_safe_test_driver

contains
subroutine linear_system_ill_defined_safe_test_driver(mpiglobal, kill_on_failure)
  !> mpi environment
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program upon failure of an assertion
  logical, intent(in), optional :: kill_on_failure

  type(unit_test_type) :: test_report
  character(len=*), parameter :: module_tested = "linear_system_ill_defined_safe"

  call test_report%init( mpiglobal )
  call test_ill_defined_safe_solve_complex_dp( test_report )
  call test_report%report( module_tested, kill_on_failure )
  call test_report%finalise()
end subroutine linear_system_ill_defined_safe_test_driver

!> Subroutine that tests [[positive_definite_solve_complex_dp]]
subroutine test_ill_defined_safe_solve_complex_dp(test_report)
  !> Test object
  type(unit_test_type), intent(inout) :: test_report

  complex(dp), parameter :: identity_2x2(2, 2) = reshape([zone, zzero, zzero, zone], [2, 2])

  call test_template( complex_positive_definite_matrix_5x5, complex_matrix_5x5, 1, test_report )
  call test_template( complex_positive_definite_matrix_5x5, complex_matrix_5x5, 2, test_report )
  call test_template( complex_positive_definite_matrix_5x5, complex_matrix_5x7, 3, test_report )
  call test_template( complex_positive_definite_matrix_5x5, complex_matrix_5x7, 4, test_report )
  call test_template( identity_2x2, reshape( complex_matrix_5x5, [2, 5] ), 5, test_report )
  call test_template( identity_2x2, reshape( complex_matrix_5x5, [2, 5] ), 6, test_report )
  contains 
    !> Template used in [[test_positive_definite_solve_complex_dp]]
    subroutine test_template( A, y, test_number, test_obj )
      !> Positive definite matrix `A` to be inverted and applied to `y`
      complex(dp), intent(in) :: A(:, :)
      !> Matrix `y`
      complex(dp), intent(in) :: y(:, :)
      !> Number to identifier this test [[positive_definite_solve]]
      integer(i32), intent(in) :: test_number
      !> Test object
      type(unit_test_type), intent(inout) :: test_obj

      real(dp), parameter :: tol = 1e-8_dp
      character(len=*), parameter :: test_identifier = "test_ill_defined_safe_solve"
      complex(dp), allocatable :: x(:, :)
      
      x = y
      ! Solve linear system
      call ill_defined_safe_solve(A, x)
      ! Check the solution: now x contains (A^-1)y
      call test_obj%assert( all_close(y, matmul(A, x), tol), &
        test_identifier // ": A*x is not equal y in test number " // to_char(test_number) // &
        ", error: " // to_char(maxval(abs(y-matmul(A, x)))) )
    end subroutine  test_template
end subroutine test_ill_defined_safe_solve_complex_dp

end module
