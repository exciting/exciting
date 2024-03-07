module gw_test_drivers
  use modmpi, only: mpiinfo

  implicit none
  
  private

  public :: gw_test_driver

contains

subroutine gw_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure 

    ! Call test drivers here

end subroutine

end module