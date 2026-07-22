!> Unit tests for ghost band detection and filtering.
module ghost_band_filter_test
  use ghost_band_filter, only: detect_ghost_bands, filter_ghost_bands, reshape_sv_arrays, restore_shape_sv_arrays
  use math_utils, only: all_close
  use modmpi, only: mpiinfo
  use precision, only: dp
  use unit_test_framework, only: unit_test_type

  implicit none

  private

  public :: run_ghost_band_filter_test_driver

contains

!> Run unit tests for ghost band filtering.
subroutine run_ghost_band_filter_test_driver(mpiglobal, kill_on_failure)
  !> MPI information.
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program before the test driver finishes if an assertion fails.
  logical, intent(in), optional :: kill_on_failure

  type(unit_test_type) :: test_report
  call test_report%init(mpiglobal)

  call test_detect_ghost_bands_across_all_columns(test_report)
  call test_filter_ghost_bands_reorders_all_columns(test_report)
  call test_filter_ghost_bands_negative_spectrum(test_report)
  call test_restore_shape_sv_arrays_selected_k_point(test_report)

  if (present(kill_on_failure)) then
    call test_report%report('ghost_band_filter', kill_on_failure)
  else
    call test_report%report('ghost_band_filter')
  end if

  call test_report%finalise()
end subroutine run_ghost_band_filter_test_driver

!> Verify that detection counts leading ghost states across all spectra.
subroutine test_detect_ghost_bands_across_all_columns(test_report)
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  real(dp), parameter :: tolerance_smallest_allowed_eval = 0.1_dp
  real(dp) :: apwe0(1, 1, 1)
  real(dp) :: lorbe0(1, 1, 1)
  real(dp) :: evals(4, 2)
  integer :: n_ghost_states

  lorbe0 = 0.0_dp
  apwe0 = 1.0_dp
  evals(:, 1) = [-0.4_dp, -0.1_dp, 0.4_dp, 0.6_dp]
  evals(:, 2) = [-0.5_dp, -0.2_dp, 0.3_dp, 0.5_dp]

  call detect_ghost_bands(lorbe0, apwe0, size(evals, dim=1), evals, &
    tolerance_smallest_allowed_eval, n_ghost_states, emit_warnings=.false.)

  call test_report%assert(n_ghost_states == 2, &
    'Ghost detection should count leading ghost bands from every passed column.')
end subroutine test_detect_ghost_bands_across_all_columns

!> Verify that filtering applies the shared shift consistently to all spectra.
subroutine test_filter_ghost_bands_reorders_all_columns(test_report)
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  complex(dp) :: evecs(1, 4, 2)
  complex(dp) :: expected_evecs(1, 4, 2)
  real(dp) :: evals(4, 2)
  real(dp) :: expected_evals(4, 2)

  evals(:, 1) = [-0.4_dp, -0.1_dp, 0.4_dp, 0.6_dp]
  evals(:, 2) = [-0.5_dp, -0.2_dp, 0.3_dp, 0.5_dp]

  evecs(1, :, 1) = cmplx([11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp], 0.0_dp, kind=dp)
  evecs(1, :, 2) = cmplx([21.0_dp, 22.0_dp, 23.0_dp, 24.0_dp], 0.0_dp, kind=dp)

  call filter_ghost_bands(size(evals, dim=1), 2, evals, evecs, emit_warnings=.false.)

  expected_evals(:, 1) = [0.4_dp, 0.6_dp, 10.6_dp, 10.6_dp]
  expected_evals(:, 2) = [0.3_dp, 0.5_dp, 10.6_dp, 10.6_dp]
  expected_evecs(1, :, 1) = cmplx([13.0_dp, 14.0_dp, 11.0_dp, 12.0_dp], 0.0_dp, kind=dp)
  expected_evecs(1, :, 2) = cmplx([23.0_dp, 24.0_dp, 21.0_dp, 22.0_dp], 0.0_dp, kind=dp)

  call test_report%assert(all_close(evals, expected_evals), &
    'Ghost filtering should move the non-ghost eigenvalues to the front in every column.')
  call test_report%assert(all_close(evecs, expected_evecs), &
    'Ghost filtering should keep eigenvectors aligned with the reordered spectra.')
end subroutine test_filter_ghost_bands_reorders_all_columns

!> Verify that shifted ghost bands still move above the spectrum when all
!> eigenvalues are negative.
subroutine test_filter_ghost_bands_negative_spectrum(test_report)
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  complex(dp) :: evecs(1, 3, 1)
  real(dp) :: evals(3, 1)
  real(dp) :: expected_evals(3, 1)

  evals(:, 1) = [-5.0_dp, -4.0_dp, -3.0_dp]
  evecs(1, :, 1) = cmplx([11.0_dp, 12.0_dp, 13.0_dp], 0.0_dp, kind=dp)

  call filter_ghost_bands(size(evals, dim=1), 1, evals, evecs, emit_warnings=.false.)

  expected_evals(:, 1) = [-4.0_dp, -3.0_dp, 47.0_dp]

  call test_report%assert(all_close(evals, expected_evals), &
    'Ghost filtering should place shifted ghost bands above the spectrum even if all eigenvalues are negative.')
end subroutine test_filter_ghost_bands_negative_spectrum

!> Verify that restoring the reshaped second-variational arrays only changes
!> the selected k-point.
subroutine test_restore_shape_sv_arrays_selected_k_point(test_report)
  !> Unit test report.
  type(unit_test_type), intent(inout) :: test_report

  complex(dp) :: evecsv(2, 2)
  complex(dp) :: evecsv_reshape(1, 2, 2)
  real(dp) :: evalsv(2, 2)
  real(dp) :: evalsv_before(2, 2)
  real(dp) :: evalsv_reshape(2, 1)

  evalsv(:, 1) = [1.0_dp, 2.0_dp]
  evalsv(:, 2) = [10.0_dp, 20.0_dp]
  evalsv_before = evalsv
  evecsv(:, 1) = cmplx([101.0_dp, 102.0_dp], 0.0_dp, kind=dp)
  evecsv(:, 2) = cmplx([201.0_dp, 202.0_dp], 0.0_dp, kind=dp)

  call reshape_sv_arrays(evalsv, evecsv, 2, 1, evalsv_reshape, evecsv_reshape)
  evalsv_reshape(:, 1) = [3.0_dp, 4.0_dp]
  evecsv_reshape(1, 1, :) = cmplx([301.0_dp, 302.0_dp], 0.0_dp, kind=dp)
  evecsv_reshape(1, 2, :) = cmplx([401.0_dp, 402.0_dp], 0.0_dp, kind=dp)

  call restore_shape_sv_arrays(evalsv_reshape, evecsv_reshape, 2, 1, evalsv, evecsv)

  call test_report%assert(all_close(evalsv(:, 1), [3.0_dp, 4.0_dp]), &
    'Restoring reshaped eigenvalues should update the selected k-point.')
  call test_report%assert(all_close(evalsv(:, 2), evalsv_before(:, 2)), &
    'Restoring reshaped eigenvalues should leave the other k-points unchanged.')
end subroutine test_restore_shape_sv_arrays_selected_k_point

end module ghost_band_filter_test
