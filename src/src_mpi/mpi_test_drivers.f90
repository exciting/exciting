module mpi_test_drivers
    use modmpi, only: mpiinfo

    use modmpi_test, only: modmpi_test_driver

    private
    public :: mpi_test_driver 

contains
    subroutine mpi_test_driver(mpiglobal, kill_on_failure)
        type(mpiinfo), intent(in) :: mpiglobal
        logical, optional :: kill_on_failure 

        call modmpi_test_driver(mpiglobal, kill_on_failure)

    end subroutine

end module mpi_test_drivers