module gw_test_drivers
  use modmpi, only: mpiinfo
  use gw_io_tests, only: run_gw_io_test_driver

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

    logical :: kill_if_fails

    kill_if_fails = .false.
    if( present(kill_on_failure) ) kill_if_fails = kill_on_failure
    
    ! Call test drivers here
    call run_gw_io_test_driver( mpiglobal, kill_if_fails )

end subroutine

end module