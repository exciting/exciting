!> Module responsible for executing all unit tests
module unit_test_drivers
   use modmpi, only: mpiinfo, terminate_mpi_env
   use cmd_line_args, only: unit_tests_type, get_unit_tests_string

   ! Load unit test driver modules here. One per src/ subdirectory 
   use math_test_drivers, only: math_test_driver
   use mpi_test_drivers, only: mpi_test_driver

   implicit none
   private
   public :: unit_test_driver

contains

   !> Run unit tests.
   !>
   !> A recommendation is to have one unit test driver per subdirectory.
   !> For example, src_gw will have one test driver, which itself calls
   !> the tests for each module within the subdirectory.
   !>
   !> If run%all = .true., all tests get run regardless of what the other
   !> logicals in run are set to.
   subroutine unit_test_driver(mpiglobal, kill_on_failure)
      !> Global MPI environment 
      type(mpiinfo), intent(inout) :: mpiglobal
      !> Immediately kill the program if an assertion fails
      logical, intent(in) :: kill_on_failure

      !> List of unit tests to run
      character(len=500) :: unit_test_list
      !> Object containing logicals, indicating which unit tests to run
      type(unit_tests_type) :: run

      unit_test_list = get_unit_tests_string(mpiglobal)

      if (len(trim(unit_test_list)) == 0) then
         call terminate_mpi_env(mpiglobal, &
            & "unit test driver has been called without passing &
            & the command line argument, -run-unit-tests")
      end if

      call run%init(unit_test_list, mpiglobal)

      if (run%math .or. run%all) then         
         call math_test_driver(mpiglobal, kill_on_failure) 
      end if

      if (run%mpi .or. run%all) then         
         call mpi_test_driver(mpiglobal, kill_on_failure)
      end if
      
      if (run%gw .or. run%all) then
         ! Placeholder
         !call gw_test_driver(mpiglobal, kill_on_failure)
      end if

   end subroutine unit_test_driver

end module unit_test_drivers
