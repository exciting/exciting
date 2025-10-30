module mpi_test_drivers
  use modmpi, only: mpiinfo
  use modmpi_test, only: modmpi_test_driver
  use exciting_mpi_test, only: exciting_mpi_test_driver

  private
  public :: mpi_test_driver 

contains
  subroutine mpi_test_driver( mpiglobal, kill_on_failure )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> Kill the program before the test driver finishes if an assertion fails
    logical, optional :: kill_on_failure 

    call modmpi_test_driver( mpiglobal, kill_on_failure )
    call exciting_mpi_test_driver( mpiglobal, kill_on_failure )
  end subroutine

end module mpi_test_drivers