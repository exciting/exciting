!> Module for collecting unit test drivers for modules in the math directory.

module math_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here
  use math_utils_test, only: math_utils_test_driver

  use grid_utils_test, only: grid_utils_test_driver

  private
  public :: math_test_driver

contains
  subroutine math_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure 

    ! Call test drivers here
    call math_utils_test_driver(mpiglobal, kill_on_failure)
    
    call grid_utils_test_driver(mpiglobal, kill_on_failure)
  end subroutine math_test_driver

end module math_test_drivers