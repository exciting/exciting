!> Module with unit tests for [[mod_core_states(module)]]
module core_states_tests
  use mod_core_states, only: build_core_state_indices
  use modmpi, only: mpiinfo
  use precision, only: i32
  use unit_test_framework, only: unit_test_type

  implicit none

  private

  public :: run_core_states_test_driver

contains

!> Run tests for [[mod_core_states(module)]]
subroutine run_core_states_test_driver( mpiglobal, kill_on_failure )
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program before the test driver finishes
  !> if an assertion fails
  logical, intent(in), optional :: kill_on_failure

  type(unit_test_type) :: test_report

  call test_report%init( mpiglobal )

  call test_contiguous_core_states( test_report )
  call test_mixed_core_states( test_report )
  call test_no_core_states( test_report )

  if (present(kill_on_failure)) then
    call test_report%report( 'core_states', kill_on_failure )
  else
    call test_report%report( 'core_states' )
  end if

  call test_report%finalise()

end subroutine

!> Unit tests for contiguous GW core-state flags.
subroutine test_contiguous_core_states( test_report )
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  integer(i32) :: core_indices(4)
  integer(i32) :: n_core

  call build_core_state_indices([.true., .true., .false., .false.], core_indices, n_core)

  call test_report%assert(n_core == 2_i32, 'Unexpected number of contiguous GW core states.')
  call test_report%assert(all(core_indices == [1_i32, 2_i32, 0_i32, 0_i32]), &
    'Contiguous GW core-state indices do not match the expected compact map.')
end subroutine

!> Unit tests for mixed GW core-state flags.
subroutine test_mixed_core_states( test_report )
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  integer(i32) :: core_indices(4)
  integer(i32) :: n_core

  call build_core_state_indices([.false., .true., .false., .true.], core_indices, n_core)

  call test_report%assert(n_core == 2_i32, 'Unexpected number of mixed GW core states.')
  call test_report%assert(all(core_indices == [2_i32, 4_i32, 0_i32, 0_i32]), &
    'Mixed GW core-state indices do not retain the species-state positions.')
end subroutine

!> Unit tests for species without GW core states.
subroutine test_no_core_states( test_report )
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  integer(i32) :: core_indices(3)
  integer(i32) :: n_core

  call build_core_state_indices([.false., .false., .false.], core_indices, n_core)

  call test_report%assert(n_core == 0_i32, 'Species without core states should report zero GW core states.')
  call test_report%assert(all(core_indices == [0_i32, 0_i32, 0_i32]), &
    'Species without core states should leave the GW core-state map empty.')
end subroutine

end module
