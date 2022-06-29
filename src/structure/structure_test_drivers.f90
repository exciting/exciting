!> Collect all unit test drivers for structure modules

module structure_test_drivers
   use modmpi, only: mpiinfo

   ! Load test drivers here
   use super_cell_utils_test, only: super_cell_utils_test_driver
   use unit_cell_utils_test, only: unit_cell_utils_test_driver

   private
   public :: structure_test_driver

contains

   subroutine structure_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure

      ! Call test drivers here
      call super_cell_utils_test_driver(mpiglobal, kill_on_failure)
      call unit_cell_utils_test_driver(mpiglobal, kill_on_failure)

   end subroutine structure_test_driver

end module structure_test_drivers

