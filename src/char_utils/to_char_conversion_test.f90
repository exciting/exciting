! Created by  on 24/03/2022.
!> Test [[to_char_conversion]]
module to_char_conversion_test
  use asserts, only: assert
  use modmpi, only: mpiinfo
  use precision, only: sp, dp
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
    !> Number of assertions
    integer, parameter :: n_assertions = 32

    call test_report%init(n_assertions, mpiglobal)

    call test_to_char_conversion(test_report)

    call test_vector_conversion(test_report)

    call test_matrix_conversion(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('to_char_conversion', kill_on_failure)
    else
      call test_report%report('to_char_conversion')
    end if
  end subroutine to_char_conversion_test_driver

  !> Test basic type conversion.
  subroutine test_to_char_conversion(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report

    ! Test logical input
    call test_report%assert(to_char(.true.) == "TRUE",   'Test to_char for logical. Expected: "TRUE".')
    call test_report%assert(to_char(.false.) == "FALSE", 'Test to_char for logical. Expected: "FALSE".')

    ! Test integer input
    call test_report%assert(to_char(1234) == "1234", 'Test to_char for integer. Expected: "1234".')
    call test_report%assert(to_char(1) == "1",       'Test to_char for integer. Expected: "1".')
    call test_report%assert(to_char(-4) == "-4",     'Test to_char for integer. Expected: "-4".')
    call test_report%assert(to_char(-0) == "0",      'Test to_char for integer. Expected: "0".')

    ! Test real input                                   
    call test_report%assert(to_char(3._sp)      ==  "3.000000E+00", 'Test to_char for real. Expected:  "3.000000E+00".')
    call test_report%assert(to_char(3.123_sp)   ==  "3.123000E+00", 'Test to_char for real. Expected:  "3.123000E+00".')
    call test_report%assert(to_char(-3._sp)     == "-3.000000E+00", 'Test to_char for real. Expected: "-3.000000E+00".')
    call test_report%assert(to_char(335.231_sp) ==  "3.352310E+02", 'Test to_char for real. Expected:  "3.352310E+02".')

    ! Test double input
    call test_report%assert(to_char(3._dp)      == "3.00000000000000E+00", &
            'Test to_char for double. Expected:  "3.00000000000000E+00".')
    call test_report%assert(to_char(3.123_dp)   ==  "3.12300000000000E+00", &
            'Test to_char for double. Expected:  "3.12300000000000E+00".')
    call test_report%assert(to_char(-3._dp)     == "-3.00000000000000E+00", &
            'Test to_char for double. Expected: "-3.00000000000000E+00".')
    call test_report%assert(to_char(335.231_dp) == "3.35231000000000E+02", &
            'Test to_char for double. Expected:  "3.35231000000000E+02".')

    ! CAUTION: These tests will fail if default_imag_identifier is changed
    ! Test real complex input
    call test_report%assert(to_char(cmplx(3._sp, 2._sp)) == "3.000000E+00+2.000000E+00i", &
            'Test to_char for complex. Expected: "3.000000E+00+2.000000E+00i".')
    call test_report%assert(to_char(cmplx(3.123_sp, -12.21342_sp)) == "3.123000E+00-1.221342E+01i", &
            'Test to_char for complex. Expected: "3.123000E+00-1.221342E+01i".')

    ! Test double complex input
    call test_report%assert(to_char(cmplx(3._dp, 2._dp, dp)) == "3.00000000000000E+00+2.00000000000000E+00i", &
            'Test to_char for complex. Expected: "3.00000000000000E+00+2.00000000000000E+00i".')
    call test_report%assert(to_char(cmplx(3.123_dp, -12.21342_dp, dp)) == "3.12300000000000E+00-1.22134200000000E+01i", &
            'Test to_char for complex. Expected: "3.12300000000000E+00-1.22134200000000E+01i".')
  end subroutine test_to_char_conversion

  !> Test vector conversion
  subroutine test_vector_conversion(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report

    call test_report%assert(to_char([.true.,.false., .false.]) == "[TRUE,FALSE,FALSE]", &
            'Test t_char for integer vector. Expected: "[TRUE,FALSE,FALSE]".')
    call test_report%assert(to_char([.TRUE., .true.]) == "[TRUE,TRUE]", &
            'Test t_char for integer vector. Expected: "[TRUE,TRUE]".')

    call test_report%assert(to_char([1, 2, 3, 5]) == "[1,2,3,5]", &
            'Test t_char for integer vector. Expected: "[1,2,3,5]".')
    call test_report%assert(to_char([1]) == "[1]", &
            'Test t_char for integer vector. Expected: "[1]".')
    
    call test_report%assert(to_char([1._sp, 0.2_sp, 30._sp, 1.2_sp]) == "[1.000000E+00,2.000000E-01,3.000000E+01,1.200000E+00]", &
            'Test t_char for real vector. Expected: "[1.000000E+00,2.000000E-01,3.000000E+01,1.200000E+00]".')
            
    call test_report%assert(to_char([1._dp, 0.2_dp, 30._dp, 1.2_dp]) == &
                    "[1.00000000000000E+00,2.00000000000000E-01,3.00000000000000E+01,1.20000000000000E+00]", &
            'Test t_char for double vector. Expected: &
                    "[1.00000000000000E+00,2.00000000000000E-01,3.00000000000000E+01,1.20000000000000E+00]".')
                    
    call test_report%assert(to_char([cmplx(3._sp, 2._sp), cmplx(3.4_sp, -0.1_sp), cmplx(0._sp, 50._sp)]) == &
                    "[3.000000E+00+2.000000E+00i,3.400000E+00-1.000000E-01i,0.000000E+00+5.000000E+01i]", &
            'Test t_char for complex vector. Expected: &
                    "[3.000000E+00+2.000000E+00i,3.400000E+00-1.000000E-01i,0.000000E+00+5.000000E+01i]".')
                    
    call test_report%assert(to_char([cmplx(3._dp, 2._dp, dp), cmplx(3.4_dp, -0.1_dp, dp), cmplx(0._dp, 50._dp, dp)]) ==&
                    "[3.00000000000000E+00+2.00000000000000E+00i,&
                            3.40000000000000E+00-1.00000000000000E-01i,&
                                    0.00000000000000E+00+5.00000000000000E+01i]", &
            'Test t_char for double complex vector. Expected: &
                    "[3.00000000000000E+00+2.00000000000000E+00i,&
                            3.40000000000000E+00-1.00000000000000E-01i,&
                                    0.00000000000000E+00+5.00000000000000E+01i]".')
  end subroutine test_vector_conversion

  !> Test matrix conversion.
  subroutine test_matrix_conversion(test_report)
    !> Unit test report
    type(unit_test_type) :: test_report
    
    call test_report%assert(to_char(reshape([2._dp, 3.4_dp, 5._dp, &
                                             0._dp, 50._dp, 3.5_dp], [3, 2])) ==&
                    "[[2.00000000000000E+00,3.40000000000000E+00,5.00000000000000E+00],&
                      [0.00000000000000E+00,5.00000000000000E+01,3.50000000000000E+00]]", &
            'Test t_char for double matrix. Expected: &
                    "[[2.00000000000000,0.000000000000000E+000],&
                      [3.500000000000000,50.0000000000000],&
                      [5.00000000000000,3.50000000000000]]".')

    call test_report%assert(to_char(transpose(reshape([cmplx(3._dp, 2._dp, dp),  cmplx(3.4_dp, -0.1_dp, dp), &
                                                       cmplx(0._dp, 50._dp, dp), cmplx(3.5_dp, 0.5_dp, dp)], &
                                                        [2, 2]))) ==&
                    "[[3.00000000000000E+00+2.00000000000000E+00i,0.00000000000000E+00+5.00000000000000E+01i],&
                      [3.40000000000000E+00-1.00000000000000E-01i,3.50000000000000E+00+5.00000000000000E-01i]]", &
            'Test t_char for double complex matrix. Expected: &
            "[[3.00000000000000E+00+2.00000000000000E+00i,0.00000000000000E+00+5.00000000000000E+01i],&
              [3.40000000000000E+00-1.00000000000000E-01i,3.50000000000000E+00+5.00000000000000E-01i]]".')
  end subroutine test_matrix_conversion

end module to_char_conversion_test