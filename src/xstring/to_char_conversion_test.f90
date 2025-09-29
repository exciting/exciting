!> Test [[to_char_conversion]]
module to_char_conversion_test
  use asserts, only: assert
  use math_utils, only: transpose_reshape
  use modmpi, only: mpiinfo
  use precision, only: dp, i32, long_int, sp
  use unit_test_framework, only: unit_test_type

  use to_char_conversion, only: to_char

  implicit none

  private
  public :: to_char_conversion_test_driver

  contains

   !> Run tests for the lattice module
  subroutine to_char_conversion_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> test object
    type(unit_test_type) :: test_report
    integer(i32), parameter :: n_assertions_test_to_char_conversion = 18
    integer(i32), parameter :: n_assertions_test_vector_conversion = 17
    integer(i32), parameter :: n_assertions_test_matrix_conversion = 4
    !> Number of assertions
    integer(i32), parameter :: n_assertions = n_assertions_test_to_char_conversion + &
                                              n_assertions_test_vector_conversion + &
                                              n_assertions_test_matrix_conversion

    call test_report%init(n_assertions, mpiglobal)

    call test_to_char_conversion(test_report)
    call test_vector_conversion(test_report)
    call test_matrix_conversion(test_report)

    call test_report%report('to_char_conversion', kill_on_failure)
  end subroutine to_char_conversion_test_driver

  !> Test basic type conversion.
  subroutine test_to_char_conversion(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report

    ! Test logical input
    call test_report%assert( to_char(.true.) == "TRUE", error_message("logical", to_char(.true.), "TRUE") )
    call test_report%assert( to_char(.false.) == "FALSE", error_message("logical", to_char(.false.), "FALSE") )

    ! Test integer input
    call test_report%assert( to_char(1234) == "1234", error_message("integer", to_char(1234), "1234") )
    call test_report%assert( to_char(1) == "1", error_message("integer", to_char(1), "1") )
    call test_report%assert( to_char(-4) == "-4", error_message("integer", to_char(-4), "-4") )
    call test_report%assert( to_char(-0) == "0", error_message("integer", to_char(-0), "0") )

    ! Test real input                                   
    call test_report%assert(to_char(3._sp)      ==  "3.000000E+00", error_message("real(sp)", to_char(3._sp),      "3.000000E+00"))
    call test_report%assert(to_char(3.123_sp)   ==  "3.123000E+00", error_message("real(sp)", to_char(3.123_sp),   "3.123000E+00"))
    call test_report%assert(to_char(-3._sp)     == "-3.000000E+00", error_message("real(sp)", to_char(-3._sp),     "-3.000000E+00"))
    call test_report%assert(to_char(335.231_sp) ==  "3.352310E+02", error_message("real(sp)", to_char(335.231_sp), "3.352310E+02"))

    ! Test double input
    call test_report%assert(to_char(3._dp)      == "3.00000000000000E+00",  error_message("real(dp)", to_char(3._dp),      "3.00000000000000E+00") )
    call test_report%assert(to_char(3.123_dp)   == "3.12300000000000E+00",  error_message("real(dp)", to_char(3.123_dp),   "3.12300000000000E+00") )
    call test_report%assert(to_char(-3._dp)     == "-3.00000000000000E+00", error_message("real(dp)", to_char(-3._dp),     "-3.00000000000000E+00") )
    call test_report%assert(to_char(335.231_dp) == "3.35231000000000E+02",  error_message("real(dp)", to_char(335.231_dp), "3.35231000000000E+02") )

    ! CAUTION: These tests will fail if default_imag_identifier is changed
    ! Test real complex input
    call test_report%assert(to_char(cmplx(3._sp, 2._sp, sp)) == "3.000000E+00+2.000000E+00i", &
      error_message("complex(sp)", to_char(cmplx(3._sp, 2._sp, sp)), "3.000000E+00+2.000000E+00i") )
    call test_report%assert(to_char(cmplx(3.123_sp, -12.21342_sp, sp)) == "3.123000E+00-1.221342E+01i", &
      error_message("complex(sp)", to_char(cmplx(3.123_sp, -12.21342_sp, sp)), "3.123000E+00-1.221342E+01i") )

    ! Test double complex input
    call test_report%assert(to_char(cmplx(3._dp, 2._dp, dp)) == "3.00000000000000E+00+2.00000000000000E+00i", &
      error_message("complex(dp)", to_char(cmplx(3._dp, 2._dp, dp)), "3.00000000000000E+00+2.00000000000000E+00i") )
    call test_report%assert(to_char(cmplx(3.123_dp, -12.21342_dp, dp)) == "3.12300000000000E+00-1.22134200000000E+01i", &
      error_message("complex(dp)", to_char(cmplx(3.123_dp, -12.21342_dp, dp)), "3.12300000000000E+00-1.22134200000000E+01i") )
  end subroutine test_to_char_conversion

  !> Test vector conversion
  subroutine test_vector_conversion(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report
    character(len=*), parameter :: expected_empty_array = "[,]"
    character(len=:), allocatable :: category

    category = "logical vector"
    call run_unit_test( test_report, category, to_char([.true.,.false., .false.]), "[TRUE,FALSE,FALSE]" )
    call run_unit_test( test_report, category, to_char([.TRUE., .true.]), "[TRUE,TRUE]" )
    call run_unit_test( test_report, category, to_char([logical:: ]), expected_empty_array )

    category = "integer(i32) vector"
    call run_unit_test( test_report, category, to_char([1, 2, 3, 5]), "[1,2,3,5]" )
    call run_unit_test( test_report, category, to_char([1]), "[1]" )
    call run_unit_test( test_report, category, to_char([integer(i32):: ]), expected_empty_array )

    category = "integer(long_int) vector"
    call run_unit_test( test_report, category, to_char([1_long_int, 4294967296_long_int]), "[1,4294967296]" )
    call run_unit_test( test_report, category, to_char([1_long_int]), "[1]" )
    call run_unit_test( test_report, category, to_char([integer(long_int):: ]), expected_empty_array )
    
    category = "real(sp) vector"
    call run_unit_test( test_report, category, to_char([1._sp, 0.2_sp, 30._sp, 1.2_sp]), "[1.000000E+00,2.000000E-01,3.000000E+01,1.200000E+00]" )
    call run_unit_test( test_report, category, to_char([real(sp):: ]), expected_empty_array )

    category = "real(dp) vector"
    call run_unit_test( test_report, category, to_char([1._dp, 0.2_dp, 30._dp, 1.2_dp]), "[1.00000000000000E+00,2.00000000000000E-01,3.00000000000000E+01,1.20000000000000E+00]" )
    call run_unit_test( test_report, category, to_char([real(dp):: ]), expected_empty_array )

    category = "complex(sp) vector"
    call run_unit_test( test_report, category, to_char([cmplx(3._sp, 2._sp, sp), cmplx(3.4_sp, -0.1_sp, sp), cmplx(0._sp, 50._sp, sp)]), &
      "[3.000000E+00+2.000000E+00i,3.400000E+00-1.000000E-01i,0.000000E+00+5.000000E+01i]" )
    call run_unit_test( test_report, category, to_char([complex(sp):: ]), expected_empty_array )

    category = "complex(dp) vector"
    call run_unit_test( test_report, category, to_char([cmplx(3._dp, 2._dp, dp), cmplx(3.4_dp, -0.1_dp, dp), cmplx(0._dp, 50._dp, dp)]), &
      "[3.00000000000000E+00+2.00000000000000E+00i,3.40000000000000E+00-1.00000000000000E-01i,0.00000000000000E+00+5.00000000000000E+01i]" )
    call run_unit_test( test_report, category, to_char([complex(dp):: ]), expected_empty_array )
  end subroutine test_vector_conversion

  !> Test matrix conversion.
  subroutine test_matrix_conversion(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report
    character(len=*), parameter :: expected_empty_matrix = "[[,]]"
    character(len=:), allocatable :: category
    
    category = "real(dp) matrix"
    call run_unit_test( test_report, category, to_char(reshape([2._dp, 3.4_dp, 5._dp, 0._dp, 50._dp, 3.5_dp], [3, 2])), &
      "[[2.00000000000000E+00,3.40000000000000E+00,5.00000000000000E+00],[0.00000000000000E+00,5.00000000000000E+01,3.50000000000000E+00]]" )
    call run_unit_test( test_report, category, to_char(reshape([real(dp):: ], [0,0])), expected_empty_matrix )

    category = "complex(dp) matrix"
    call run_unit_test( test_report, category, to_char( &
      transpose_reshape([cmplx(3._dp, 2._dp, dp),  cmplx(3.4_dp, -0.1_dp, dp), &
                         cmplx(0._dp, 50._dp, dp), cmplx(3.5_dp, 0.5_dp, dp)], [2, 2]) ), &
                   &"[[3.00000000000000E+00+2.00000000000000E+00i,0.00000000000000E+00+5.00000000000000E+01i],&
                     &[3.40000000000000E+00-1.00000000000000E-01i,3.50000000000000E+00+5.00000000000000E-01i]]" )
    call run_unit_test( test_report, category, to_char(reshape([complex(dp):: ], [0,0])), expected_empty_matrix )
  end subroutine test_matrix_conversion

  !> (private) Run a unit test, checking if an observed result corresponds to the expected one
  subroutine run_unit_test( test_report, category, observed, expected ) 
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report
    !> String with the category (kind) where the test failed
    character(len=*), intent(in) :: category
    !> Observed result
    character(len=*), intent(in) :: observed
    !> Expected result
    character(len=*), intent(in) :: expected
    
    call test_report%assert( observed == expected, error_message(category, observed, expected) )
  end subroutine

  !> (private) Returns the error message for a failing test
  pure function error_message( category, observed, expected ) result(msg)
    !> String with the category (kind) where the test failed
    character(len=*), intent(in) :: category
    !> Observed result
    character(len=*), intent(in) :: observed
    !> Expected result
    character(len=*), intent(in) :: expected
    !> Resulting message
    character(len=:), allocatable :: msg
    msg = "Test to_char for " // category // " failed. Result: " // observed // ". Expected: " // expected
  end function

end module to_char_conversion_test
