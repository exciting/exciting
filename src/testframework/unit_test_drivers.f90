!> Module responsible for executing all unit tests
module unit_test_drivers
   use modmpi, only: mpiinfo, terminate_mpi_env
   use cmd_line_args, only: unit_tests_type, get_unit_tests_string

   ! Load unit test driver modules here

   implicit none
   private
   public :: unit_test_driver

contains

   !> Setup for unit tests.
   subroutine setup()
      use modmpi, only: initmpi
      call initmpi
   end subroutine

   !> Tear down for unit tests.
   subroutine teardown()
      use modmpi, only: finitmpi
      call finitmpi
   end subroutine

   !> Run unit tests.
   !>
   !> A recommendation is to have one unit test driver per subdirectory.
   !> For example, src_gw will have one test driver, which itself calls
   !> the tests for each module within the subdirectory.
   !>
   !> If run%all = .true., all tests get run regardless of what the other
   !> logicals in run are set to.
   !>
   !> If one is writing tests using MPI, and they need to initialise/finalise
   !> mpi_comm_world, this routine can be modified such that setup and teardown
   !> are not called in those instances. In which case, one may need to make
   !> mpiglobal an optional arg for run%init, and if it is not passed,
   !> get mpi_comm_world from the mpi module.
   subroutine unit_test_driver(kill_on_failure)
      use modmpi, only: mpiglobal
      !> Immediately kill the program upon a failed assertion
      logical, intent(in) :: kill_on_failure

      !> List of unit tests to run
      character(len=500) :: unit_test_list
      !> Unit tests to run
      type(unit_tests_type) :: run

      ! Initialise mpiglobal
      call setup()

      unit_test_list = get_unit_tests_string()

      if (len(trim(unit_test_list)) == 0) then
         call terminate_mpi_env(mpiglobal, &
            & "unit test driver has been called without passing &
            & the command line argument, -run-unit-tests")
      end if

      call run%init(unit_test_list, mpiglobal)

      if (run%math_utils .or. run%all) then
         ! Placeholder
         !call math_utils_test_driver(mpiglobal, kill_on_failure)
      end if
      if (run%gw .or. run%all) then
         ! Placeholder
         !call gw_test_driver(mpiglobal, kill_on_failure)
      end if

      ! Finalise mpiglobal
      call teardown()

   end subroutine unit_test_driver

end module unit_test_drivers
