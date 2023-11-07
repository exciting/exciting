!> Module responsible for executing all unit tests
module unit_test_drivers
   use modmpi, only: mpiinfo, terminate_mpi_env
   use cmd_line_args, only: get_unit_tests_string
   use unit_tests, only: unit_tests_type

   ! Load unit test driver modules here. One per src/ subdirectory 
   use math_test_drivers, only: math_test_driver
   use xs_test_drivers, only: xs_test_driver
   use lapack_wrappers_test_drivers, only: lapack_wrappers_test_driver
   use mpi_test_drivers, only: mpi_test_driver
   use structure_test_drivers, only: structure_test_driver
   use testframework_test_drivers, only: testframework_test_driver
   use xstring_test_drivers, only: xstring_test_driver
   use file_io_test_drivers, only: file_io_test_driver
   use simplified_input_test_drivers, only: simplified_input_test_driver
   use hybrids_test_drivers, only: hybrids_test_driver
   use matrix_elements_test_drivers, only: matrix_elements_test_driver
   use xgrid_test_drivers, only: xgrid_test_driver
   use xhdf5_test, only: xhdf5_test_driver

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

      if (run%xs .or. run%all) then
         call xs_test_driver(mpiglobal, kill_on_failure) 
      end if
         
      if (run%xs .or. run%all) then
         call xs_test_driver(mpiglobal, kill_on_failure) 
      end if

      if (run%lapack .or. run%all) then
         call lapack_wrappers_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%mpi .or. run%all) then
         call mpi_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%structure .or. run%all) then
         call structure_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%gw .or. run%all) then
         ! Placeholder
         !call gw_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%xstring .or. run%all) then
         call xstring_test_driver(mpiglobal, kill_on_failure) 
      end if
      
      if (run%simplified_input .or. run%all) then         
         call simplified_input_test_driver(mpiglobal, kill_on_failure)
      end if 
      
      if (run%file_io .or. run%all) then
        call file_io_test_driver(mpiglobal, kill_on_failure)
      end if
      
      if (run%matrix_elements .or. run%all) then
         call matrix_elements_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%testframework .or. run%all) then
         call testframework_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%hybrids .or. run%all) then         
         call hybrids_test_driver(mpiglobal, kill_on_failure) 
      end if

      if (run%xgrid .or. run%all) then
         call xgrid_test_driver(mpiglobal, kill_on_failure)
      end if

      if (run%xhdf5 .or. run%all) then
         call xhdf5_test_driver(mpiglobal, kill_on_failure)
      end if
   end subroutine unit_test_driver

end module unit_test_drivers
