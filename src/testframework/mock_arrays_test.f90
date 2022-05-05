module mock_arrays_test

    use unit_test_framework, only: unit_test_type
    use precision, only: sp, dp
    use math_utils, only: all_zero, all_close
    use modmpi, only: mpiinfo
    use mock_arrays, only: fill_array, value_map_real_rank1,&
                        value_map_complex_rank2, &
                        value_map_complex_rank3, &
                        value_map_complex_rank4


    implicit none

    private

    public :: mock_arrays_test_driver

contains

    !> Runs all tests in this module
    !> Called from top-level routine
    subroutine mock_arrays_test_driver(mpiglobal, kill_on_failure)

        !> mpi information
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program before the test driver finishes,
        !> if an assertion fails
        logical, optional :: kill_on_failure

        !> Our test object that looks like the Zofu object
        type(unit_test_type) :: test_report
        !> Number of tests
        integer, parameter :: n_assertions = 4

        ! Initialise tests
        call test_report%init(n_assertions, mpiglobal)

        ! Call tests

        call test_fill_real_array_rank1(test_report)

        call test_fill_complex_array_rank2(test_report)

        call test_fill_complex_array_rank3(test_report)

        call test_fill_complex_array_rank4(test_report)

        ! report results
        if (present(kill_on_failure)) then
            call test_report%report('mock_arrays', kill_on_failure)
        else
            call test_report%report('mock_arrays')
        end if

        ! Deallocate after tests
        call test_report%finalise()
    end subroutine mock_arrays_test_driver

    
    !> Test fill_array for real array of rank 1.
    subroutine test_fill_real_array_rank1(test_report)
        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Real test input
        real(dp) :: a(5)

        ! Fill array with (arbitray) numbers
        call fill_array(a, value_map_real_rank1)

        call test_report%assert(all_close(a(5), 7.5_dp), &
                                'Test fill_array for real array of rank 1: &
                                Expected: Element given by a_5 = 7.5.')

    end subroutine test_fill_real_array_rank1


    !> Test fill_array for complex array of rank 2.
    subroutine test_fill_complex_array_rank2(test_report)
        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Real test input
        complex(dp) :: a(5, 5)

        ! Fill array with (arbitray) numbers
        call fill_array(a, value_map_complex_rank2)
        
        call test_report%assert(all_close(a(3, 5), cmplx(3, 2*5, dp)), &
                                'Test fill_array for complex array of rank 2: &
                                Expected: Element given by a_{3,5} = 3 + 10*i.')

    end subroutine test_fill_complex_array_rank2



    !> Test fill_array for complex array of rank 3.
    subroutine test_fill_complex_array_rank3(test_report)
        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Real test input
        complex(dp) :: a(5, 5, 5)

        ! Fill array with numbers according to value map
        call fill_array(a, value_map_complex_rank3)
        
        call test_report%assert(all_close(a(3, 4, 5), cmplx(3, 2*(4+5), dp)), &
                                'Test fill_array for complex array of rank 3: &
                            Expected: Element given by a_{3,6,7} = 3 + 20*i.')

    end subroutine test_fill_complex_array_rank3


    !> Test fill_array for complex array of rank 4.
    subroutine test_fill_complex_array_rank4(test_report)
        !> Unit test report
        type(unit_test_type), intent(inout) :: test_report
        !> tolerance
        real(dp), parameter :: tol = 1.e-10_dp
        !> Real test input
        complex(dp) :: a(3, 3, 3, 4)

        ! Fill array with numbers according to value map
        call fill_array(a, value_map_complex_rank4)
        
        call test_report%assert(all_close(a(1, 2, 3, 4), cmplx(1+3, 2*2*4, dp)), &
                                'Test fill_array for complex array of rank 4: &
                                Expected: Element given by a_{1,2,3,4} = 4 + 16*i.')

    end subroutine test_fill_complex_array_rank4

end module mock_arrays_test
