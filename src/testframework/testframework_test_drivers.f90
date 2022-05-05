!> Module for collecting unit test drivers for modules in the 
!> testframework directory. 
module testframework_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here
  use mock_arrays_test, only: mock_arrays_test_driver

  private
  public :: testframework_test_driver

contains
  subroutine testframework_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure 

    ! Call test drivers here
    call mock_arrays_test_driver(mpiglobal, kill_on_failure)

  end subroutine testframework_test_driver

end module testframework_test_drivers