!> Module with tests for the [[add_scaled_array]] module
module add_scaled_array_test
  use add_scaled_array, only: scaled_add
  use constants, only: real_one, real_zero, zi, zone, zzero
  use exciting_mpi, only: mpiinfo
  use math_utils, only: all_close
  use mock_arrays, only: complex_vector_5, complex_vector_7, complex_matrix_5x7, &
    complex_matrix_7x5, complex_matrix_5x5
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none

  private

  public :: add_scaled_array_test_driver

contains
  subroutine add_scaled_array_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional, intent(in) :: kill_on_failure
    
    type(unit_test_type) :: test_report
    integer(i32), parameter :: n_assertions_test_scaled_add_rank1_arrays = 3
    integer(i32), parameter :: n_assertions_test_scaled_add_rank2_arrays = 3
    integer(i32), parameter :: n_assertions_test_scaled_add_rank3_arrays = 6
    integer(i32), parameter :: n_assertions = n_assertions_test_scaled_add_rank1_arrays + &
                                              n_assertions_test_scaled_add_rank2_arrays + &
                                              n_assertions_test_scaled_add_rank3_arrays

    character(len=*), parameter :: module_tested = "add_scaled_array"

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_scaled_add_rank1_arrays( test_report )
    call test_scaled_add_rank2_arrays( test_report )
    call test_scaled_add_rank3_arrays( test_report )

    call test_report%report( module_tested, kill_on_failure )

    call test_report%finalise()

  end subroutine

  subroutine test_scaled_add_rank1_arrays( test_report )
    !> The test object
    type(unit_test_type), intent(inout) :: test_report

    integer(i32) :: test_number
    character(len=*), parameter :: test_id = "test_scaled_add_rank1_arrays"
    complex(dp) :: a
    complex(dp), allocatable :: x(:), y(:), y_expected(:)
    
    test_number = 1
    a = cmplx(1.0_dp, 2.0_dp, dp)
    x = complex_vector_5
    y = complex_vector_7(1:size(x))
    y_expected = a*x + y
    call scaled_add( a, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    call scaled_add( zzero, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    y = complex_vector_7(1:size(x))
    y_expected = x + y
    call scaled_add( zone, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )
  end subroutine

  subroutine test_scaled_add_rank2_arrays( test_report )
    !> The test object
    type(unit_test_type), intent(inout) :: test_report

    integer(i32) :: test_number
    character(len=*), parameter :: test_id = "test_scaled_add_rank2_arrays"
    complex(dp) :: a
    complex(dp), allocatable :: x(:, :), y(:, :), y_expected(:, :)
    
    test_number = 1
    a = cmplx(1.0_dp, 2.333_dp, dp)
    x = complex_matrix_5x5
    y = reshape( complex_matrix_7x5, shape(x) )
    y_expected = a*x + y
    call scaled_add( a, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    call scaled_add( zzero, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    x = complex_matrix_5x7
    y = zi + transpose( complex_matrix_7x5 )
    y_expected = x + y
    call scaled_add( zone, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )
  end subroutine

  subroutine test_scaled_add_rank3_arrays( test_report )
    !> The test object
    type(unit_test_type), intent(inout) :: test_report

    integer(i32) :: test_number
    character(len=*), parameter :: test_id = "test_scaled_add_rank3_arrays"
    real(dp) :: alpha
    complex(dp) :: a
    complex(dp), parameter :: M(3, 2, 5) = reshape( complex_matrix_7x5, shape(M) )
    complex(dp), parameter :: N(3, 2, 5) = reshape( zi*complex_matrix_5x7, shape(N) )
    complex(dp), allocatable :: x(:, :, :), y(:, :, :), y_expected(:, :, :)
    
    test_number = 1
    a = cmplx(1.0_dp, 2.333_dp, dp)
    x = M
    y = N
    y_expected = a*x + y
    call scaled_add( a, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    call scaled_add( zzero, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    x = M + (5._dp + zi)**N
    y = zone + M
    y_expected = x + y
    call scaled_add( zone, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    ! Tests with a real scaling factor
    test_number = test_number + 1
    alpha = 5.9873412_dp
    x = M
    y = N
    y_expected = alpha*x + y
    call scaled_add( alpha, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    call scaled_add( real_zero, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )

    test_number = test_number + 1
    x = M + (5._dp + zi)**N
    y = zone + M
    y_expected = x + y
    call scaled_add( real_one, x, y )
    call test_report%assert( all_close( y, y_expected ), report_message( test_id, "complex", test_number ) )
  end subroutine

  !> (private) Generate a report message, given a test case and test number
  function report_message( test_id, test_case, test_number ) result(message)
    !> Name of the subroutine calling `report_message`
    character(len=*), intent(in) :: test_id
    !> Test (sub)identifier, pointing out which case in [[test_id]] is tested
    character(len=*), intent(in) :: test_case
    !> Number for the test
    integer(i32), intent(in) :: test_number
    !> Message to return
    character(len=:), allocatable :: message
    message = test_id // ' - ' // test_case // ' does not match reference: test - ' // to_char(test_number)
  end function
  
end module