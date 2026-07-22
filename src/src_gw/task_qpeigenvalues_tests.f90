!> Unit tests for [[task_QPEigenvalues(module)]] helpers.
module task_qpeigenvalues_tests
  use modmpi, only: mpiinfo
  use precision, only: dp, i32
  use task_QPEigenvalues, only: index_CBm, index_VBM
  use unit_test_framework, only: unit_test_type

  implicit none

  private

  public :: run_task_qpeigenvalues_test_driver

contains

  !> Run tests for `task_QPEigenvalues`.
  subroutine run_task_qpeigenvalues_test_driver(mpiglobal, kill_on_failure)
    !> MPI information.
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion.
    logical, intent(in), optional :: kill_on_failure

    type(unit_test_type) :: test_report

    call test_report%init(mpiglobal)

    call test_band_edge_indexes_use_global_extrema(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('task_QPEigenvalues', kill_on_failure)
    else
      call test_report%report('task_QPEigenvalues')
    end if

    call test_report%finalise()
  end subroutine run_task_qpeigenvalues_test_driver

  !> Test that band-edge indexes come from the global band-edge energy.
  subroutine test_band_edge_indexes_use_global_extrema(test_report)
    !> Unit test report.
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: conduction_eigenvalues(4, 2)
    real(dp) :: valence_eigenvalues(4, 2)
    real(dp), parameter :: e_fermi = 0.0_dp

    valence_eigenvalues(:, 1) = [-4.0_dp, -1.0_dp, 1.0_dp, 2.0_dp]
    valence_eigenvalues(:, 2) = [-4.0_dp, -3.0_dp, -2.0_dp, 1.0_dp]

    conduction_eigenvalues(:, 1) = [-2.0_dp, 1.0_dp, 3.0_dp, 4.0_dp]
    conduction_eigenvalues(:, 2) = [-2.0_dp, -1.0_dp, 0.1_dp, 4.0_dp]

    call test_report%assert(index_VBM(valence_eigenvalues, e_fermi) == 2_i32, &
      'VBM index should identify the band of the global valence-band maximum.')
    call test_report%assert(index_CBm(conduction_eigenvalues, e_fermi) == 3_i32, &
      'CBM index should identify the band of the global conduction-band minimum.')
  end subroutine test_band_edge_indexes_use_global_extrema

end module task_qpeigenvalues_tests
