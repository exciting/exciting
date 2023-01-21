!> Module for collecting unit test drivers for modules in the xs directory.
module xs_test_drivers
    use modmpi, only: mpiinfo
    ! Load test drivers here
    use phonon_screening_tests, only: phonon_screening_test_driver
    use putgeteps0_tests, only: putgeteps0_test_driver
    use expand_add_eps_tests, only: expand_add_eps_test_driver

    private
    public :: xs_test_driver

contains
    subroutine xs_test_driver(mpiglobal, kill_on_failure)
        !> mpi information
        type(mpiinfo), intent(in) :: mpiglobal
        !> Kill the program before the test driver finishes
        !> if an assertion fails
        logical, optional :: kill_on_failure

        call phonon_screening_test_driver(mpiglobal, kill_on_failure)

        if (mpiglobal%is_root) then
            call putgeteps0_test_driver(mpiglobal, kill_on_failure)
        end if

        call expand_add_eps_test_driver(mpiglobal, kill_on_failure)

        
    end subroutine xs_test_driver

end module xs_test_drivers
