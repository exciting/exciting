!> Module for unit tests for the functions in [[os_utils(module)]].
module os_utils_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use os_utils

  implicit none
  private

  character(64), parameter :: TEST_DIR = 'OS_UTILS_TEST_DIR/'

  public :: run_os_utils_test_driver

  contains

    !> Run tests for [[os_utils(module)]]
    subroutine run_os_utils_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(inout) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure
      !> Test report object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 8

      character(:), allocatable :: c

      ! Initialize test object
      call test_report%init(n_assertions, mpiglobal)

      ! Run and assert tests
      call test_make_dir(test_report, mpiglobal)
      call test_remove_dir(test_report, mpiglobal)
      call test_join_paths(test_report)

      ! report results
      if (present(kill_on_failure)) then
        call test_report%report('os_utils', kill_on_failure)
      else
        call test_report%report('os_utils')
      end if

      ! Finalise test object
      call test_report%finalise()
    end subroutine run_os_utils_test_driver

    !> Test `make_directory`
    subroutine test_make_dir(test_report, comm)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(inout) :: comm

      integer :: ierr
      logical :: exists

      ierr = make_directory(TEST_DIR, comm)

      call test_report%assert( ierr == 0, &
        'Expected: Return value 0 in `make_directory`.')

      call test_report%assert( path_exists(TEST_DIR, ierr), &
        'Expected: Test directory does exist.')
    end subroutine test_make_dir

    !> Test `remove_directory`
    subroutine test_remove_dir(test_report, comm)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report
      !> MPI communicator
      type(mpiinfo), intent(inout) :: comm

      integer :: ierr

      ierr = remove_directory(TEST_DIR, comm)

      call test_report%assert( ierr == 0, &
        'Expected: Return value 0 in `remove_directory`.')

      call test_report%assert( .not. path_exists(TEST_DIR, ierr), &
        'Expected: Test directory does not exist.')
    end subroutine test_remove_dir

    !> Test 'join paths'
    subroutine test_join_paths(test_report)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report

      call test_report%assert(join_paths('path/one/', 'path/two') == 'path/one/path/two', &
              "join_paths('path/one/', 'path/two') /= 'path/one/path/two'")

      call test_report%assert(join_paths('path/one', 'path/two') == 'path/one/path/two', &
              "join_paths('path/one', 'path/two') /= 'path/one/path/two'")

      call test_report%assert(join_paths('path/one', '/path/two') == 'path/one/path/two', &
              "join_paths('path/one', '/path/two') /= 'path/one/path/two'")

      call test_report%assert(join_paths('path/one/', '/path/two') == 'path/one/path/two', &
              "join_paths('path/one/', '/path/two') /= 'path/one/path/two'")

      ! Todo(Bene #161): Fix this test.
      !call test_report%assert(join_paths('path/one////', '///path/two') == 'path/one/path/two', &
      !        "join_paths('path/one////', '///path/two') /= 'path/one/path/two'")
    end subroutine test_join_paths

end module os_utils_test