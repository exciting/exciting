!> Unit tests for [[evgw0_validation(module)]].
module evgw0_validation_tests
  use evgw0_validation, only: has_checkpoint_qp_energy_for_band, evgw0_seed_eqpsolver_is_supported, &
    should_warn_about_ignored_evgw0_eqpsolver, task_group_restart_files_exist
  use file_utils, only: delete_file
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use mod_kqpts, only: has_full_k_point_coverage
  use modmpi, only: mpiinfo
  use os_utils, only: join_paths, system_cmd
  use precision, only: i32, dp, long_int
  use quasiparticle_energies, only: checkpoint_matches_expected_qp_window
  use to_char_conversion, only: to_char
  use unit_test_framework, only: unit_test_type

  implicit none

  private

  public :: run_evgw0_validation_test_driver

contains

  !> Run unit tests for the explicit evGW0 validation helpers.
  subroutine run_evgw0_validation_test_driver(mpiglobal, kill_on_failure)
    !> MPI information.
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion.
    logical, intent(in), optional :: kill_on_failure

    type(unit_test_type) :: test_report

    call test_report%init(mpiglobal)

    call test_has_full_k_point_coverage(test_report)
    call test_has_checkpoint_qp_energy_for_band(test_report)
    call test_evgw0_eqpsolver_policy(test_report)
    call test_checkpoint_matches_expected_qp_window(test_report, mpiglobal%rank)
    call test_task_group_restart_files_exist(test_report, mpiglobal%rank)

    call test_report%report('evgw0_validation', kill_on_failure)

    call test_report%finalise()
  end subroutine run_evgw0_validation_test_driver

  !> Test the full-k-point coverage helper.
  subroutine test_has_full_k_point_coverage(test_report)
    !> Unit test report.
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(has_full_k_point_coverage([1_i32, 2_i32], 2_i32), &
      'Full irreducible-k-point coverage should be detected when all k-points are selected.')
    call test_report%assert(.not. has_full_k_point_coverage([1_i32], 2_i32), &
      'Partial k-point coverage must not be accepted for explicit evGW0.')
    call test_report%assert(.not. has_full_k_point_coverage([1_i32, 1_i32], 2_i32), &
      'Duplicate k-points must not be accepted as full explicit evGW0 coverage.')
    call test_report%assert(.not. has_full_k_point_coverage([1_i32, 3_i32], 2_i32), &
      'Out-of-range k-points must not be accepted as full explicit evGW0 coverage.')
  end subroutine test_has_full_k_point_coverage

  !> Test the helper that decides whether a band is covered by a stored evGW0 checkpoint.
  subroutine test_has_checkpoint_qp_energy_for_band(test_report)
    !> Unit test report.
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(has_checkpoint_qp_energy_for_band(4_i32, 1_i32, 6_i32), &
      'Bands inside the stored evGW0 checkpoint window should reuse checkpoint quasiparticle energies.')
    call test_report%assert(.not. has_checkpoint_qp_energy_for_band(7_i32, 1_i32, 6_i32), &
      'Bands above the stored evGW0 checkpoint window must fall back to KS energies.')
    call test_report%assert(.not. has_checkpoint_qp_energy_for_band(0_i32, 1_i32, 6_i32), &
      'Bands below the stored evGW0 checkpoint window must fall back to KS energies.')
  end subroutine test_has_checkpoint_qp_energy_for_band

  !> Test the explicit evGW0 policy for seed solvers and reused-iteration warnings.
  subroutine test_evgw0_eqpsolver_policy(test_report)
    !> Unit test report.
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(evgw0_seed_eqpsolver_is_supported(0_i32), &
      'Explicit evGW0 should accept eqpsolver=0 for the initial full iteration.')
    call test_report%assert(evgw0_seed_eqpsolver_is_supported(1_i32), &
      'Explicit evGW0 should accept eqpsolver=1 for the initial full iteration.')
    call test_report%assert(evgw0_seed_eqpsolver_is_supported(2_i32), &
      'Explicit evGW0 should accept eqpsolver=2 for the initial full iteration.')
    call test_report%assert(.not. evgw0_seed_eqpsolver_is_supported(3_i32), &
      'Explicit evGW0 should reject unsupported eqpsolver values.')
    call test_report%assert(.not. should_warn_about_ignored_evgw0_eqpsolver(0_i32, .true.), &
      'Explicit evGW0 should not warn about ignored eqpsolver when eqpsolver=0.')
    call test_report%assert(should_warn_about_ignored_evgw0_eqpsolver(1_i32, .true.), &
      'Explicit evGW0 should warn when reused iterations ignore eqpsolver=1.')
    call test_report%assert(should_warn_about_ignored_evgw0_eqpsolver(2_i32, .true.), &
      'Explicit evGW0 should warn when reused iterations ignore eqpsolver=2.')
    call test_report%assert(.not. should_warn_about_ignored_evgw0_eqpsolver(2_i32, .false.), &
      'Explicit evGW0 should not warn before checkpoint quasiparticle energies are reused.')
  end subroutine test_evgw0_eqpsolver_policy

  !> Test the evGW0 checkpoint metadata validation used for resumed runs.
  subroutine test_checkpoint_matches_expected_qp_window(test_report, mpi_rank)
    !> Unit test report.
    type(unit_test_type), intent(inout) :: test_report
    !> MPI rank running the test.
    integer(i32), intent(in) :: mpi_rank

    character(len=:), allocatable :: checkpoint_file
    real(dp) :: kpoint_vectors(3, 2)
    real(dp) :: mismatched_kpoint_vectors(3, 2)

    checkpoint_file = '/tmp/evgw0_validation_checkpoint_rank' // to_char(mpi_rank) // '.OUT'
    kpoint_vectors(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    kpoint_vectors(:, 2) = [0.25_dp, 0.0_dp, 0.0_dp]
    mismatched_kpoint_vectors = kpoint_vectors
    mismatched_kpoint_vectors(:, 2) = [0.5_dp, 0.0_dp, 0.0_dp]

    call write_mock_checkpoint(checkpoint_file, kpoint_vectors, 1_i32, 6_i32)

    call test_report%assert(checkpoint_matches_expected_qp_window(checkpoint_file, kpoint_vectors, 1_i32, 6_i32), &
      'Explicit evGW0 restart should accept a checkpoint with matching k-point and band metadata.')
    call test_report%assert(.not. checkpoint_matches_expected_qp_window(checkpoint_file, mismatched_kpoint_vectors, 1_i32, 6_i32), &
      'Explicit evGW0 restart should reject a checkpoint with mismatched k-point vectors.')
    call test_report%assert(.not. checkpoint_matches_expected_qp_window(checkpoint_file, kpoint_vectors, 1_i32, 7_i32), &
      'Explicit evGW0 restart should reject a checkpoint with mismatched upper band metadata.')

    call remove_if_present(checkpoint_file)
  end subroutine test_checkpoint_matches_expected_qp_window

  !> Test the boolean restart prerequisite probe for split `taskGroup` evGW0 jobs.
  subroutine test_task_group_restart_files_exist(test_report, mpi_rank)
    !> Unit test report.
    type(unit_test_type), intent(inout) :: test_report
    !> MPI rank running the test.
    integer(i32), intent(in) :: mpi_rank

    real(dp) :: qpoint_vectors(3, 1)
    character(len=40), parameter :: file_names(8) = [character(len=40) :: &
      'SGI_Q1.OUT', 'BARC_Q1.OUT', 'INVERSE-EPSILON_Q1.OUT', 'INVERSE-EPSH.OUT', &
      'INVERSE-EPSW1.OUT', 'INVERSE-EPSW2.OUT', 'VXCNN.OUT', 'SIGMAX_K2.OUT']
    integer(i32) :: ierr
    integer(i32) :: ifile
    integer(i32) :: file_unit
    character(len=:), allocatable :: test_directory
    character(len=:), allocatable :: file_path

    qpoint_vectors(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    test_directory = '/tmp/evgw0_task_group_restart_probe_rank' // to_char(mpi_rank)

    ierr = system_cmd('rm -rf ' // trim(test_directory))
    call test_report%assert(ierr == 0, &
      'Split evGW0 task-group restart probe test should be able to remove its temporary directory.')
    ierr = system_cmd('mkdir -p ' // trim(test_directory))
    call test_report%assert(ierr == 0, &
      'Split evGW0 task-group restart probe test should be able to create its temporary directory.')

    call test_report%assert(.not. task_group_restart_files_exist(output_format='binary', is_task_sigmac=.true., &
      is_task_qpeigenvalues=.true., qpoint_vectors=qpoint_vectors, qpeigenvalues_kpoint_indexes=[2_i32], &
      directory=test_directory), &
      'Split evGW0 task-group restart probe should reject an incomplete file set.')

    do ifile = 1, size(file_names)
      file_path = join_paths(test_directory, trim(file_names(ifile)))
      open(newunit=file_unit, file=trim(file_path), action='write', status='replace')
      close(file_unit)
    end do

    call test_report%assert(task_group_restart_files_exist(output_format='binary', is_task_sigmac=.true., &
      is_task_qpeigenvalues=.true., qpoint_vectors=qpoint_vectors, qpeigenvalues_kpoint_indexes=[2_i32], &
      directory=test_directory), &
      'Split evGW0 task-group restart probe should accept a complete file set.')

    ierr = system_cmd('rm -rf ' // trim(test_directory))
    call test_report%assert(ierr == 0, &
      'Split evGW0 task-group restart probe test should clean up its temporary directory.')
  end subroutine test_task_group_restart_files_exist

  !> Write a mock `EVALQP` checkpoint compatible with the evGW0 metadata probe.
  subroutine write_mock_checkpoint(file_name, kpoint_vectors, first_qp_band, last_qp_band)
    !> Output file name.
    character(len=*), intent(in) :: file_name
    !> Stored irreducible k-point vectors in lattice coordinates.
    real(dp), intent(in) :: kpoint_vectors(:, :)
    !> Stored first QP band.
    integer(i32), intent(in) :: first_qp_band
    !> Stored last QP band.
    integer(i32), intent(in) :: last_qp_band

    real(dp), allocatable :: eqp(:), eks(:)
    real(dp) :: efks, efqp
    integer(i32) :: ik, nkpt, unit_number
    integer(long_int) :: record_length

    call remove_if_present(file_name)
    nkpt = size(kpoint_vectors, 2)
    efks = 0.0_dp
    efqp = 0.0_dp
    allocate(eqp(first_qp_band:last_qp_band), eks(first_qp_band:last_qp_band))
    eqp = 0.0_dp
    eks = 0.0_dp

    call inquire_large(record_length, [nkpt, first_qp_band, last_qp_band], kpoint_vectors(:, 1), eqp, eks, [efqp, efks])
    call open_direct_unformatted_large(unit_number, trim(file_name), "write", record_length, "replace")
    do ik = 1, nkpt
      eqp = real(ik, dp)
      eks = -real(ik, dp)
      write(unit_number, rec=ik) nkpt, first_qp_band, last_qp_band, kpoint_vectors(:, ik), eqp, eks, efqp, efks
    end do
    close(unit_number)
    deallocate(eqp, eks)
  end subroutine write_mock_checkpoint

  !> Remove a file if it exists.
  subroutine remove_if_present(file_name)
    !> File to remove.
    character(len=*), intent(in) :: file_name

    integer(i32) :: ierr

    call delete_file(file_name, ierr)
  end subroutine remove_if_present

end module evgw0_validation_tests
