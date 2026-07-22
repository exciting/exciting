module gw_test_drivers
  use core_states_tests, only: run_core_states_test_driver
  use modmpi, only: mpiinfo
  use evgw0_validation_tests, only: run_evgw0_validation_test_driver
  use gw_io_tests, only: run_gw_io_test_driver
  use task_qpeigenvalues_tests, only: run_task_qpeigenvalues_test_driver

  implicit none

  private

  public :: gw_test_driver

contains

!> Call all unit tests related to GW.
subroutine gw_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, intent(in), optional :: kill_on_failure

    logical :: kill_if_fails

    kill_if_fails = .false.
    if( present(kill_on_failure) ) kill_if_fails = kill_on_failure

    ! Call test drivers here
    call run_core_states_test_driver( mpiglobal, kill_if_fails )
    call run_gw_io_test_driver( mpiglobal, kill_if_fails )
    call run_evgw0_validation_test_driver( mpiglobal, kill_if_fails )
    call run_task_qpeigenvalues_test_driver( mpiglobal, kill_if_fails )

end subroutine

end module
