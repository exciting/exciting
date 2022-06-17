! Created by  on 24/03/2022.

module char_utils_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here
  use to_char_conversion_test, only: to_char_conversion_test_driver

  implicit none

  private
  public :: char_utils_test_driver

  contains

  subroutine char_utils_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure

    call to_char_conversion_test_driver(mpiglobal, kill_on_failure)
  end subroutine char_utils_test_driver

end module char_utils_test_drivers