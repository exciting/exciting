!> Helper routines for validating the explicit evGW0 workflow.
module evgw0_validation
#include "asserts.fpp"
  use gw_io, only: build_file_name
  use mod_misc_gw, only: gammapoint
  use mod_selfenergy, only: file_name_sigmac, file_name_sigmax
  use mod_vxc, only: get_vxcnn_file_name
  use modmpi, only: terminate_if_false
  use os_utils, only: join_paths
  use precision, only: i32, dp
  use quasiparticle_energies, only: checkpoint_matches_expected_qp_window

  implicit none

  private

  public :: assert_checkpoint_matches_expected_qp_window
  public :: assert_task_group_restart_files_exist
  public :: evgw0_seed_eqpsolver_is_supported
  public :: has_checkpoint_qp_energy_for_band
  public :: should_warn_about_ignored_evgw0_eqpsolver
  public :: task_group_restart_files_exist

contains

  !> Return whether an evGW0 checkpoint contains a quasiparticle energy for the given band.
  pure logical function has_checkpoint_qp_energy_for_band(band_index, first_qp_band, last_qp_band) result(flag)
    !> Band index requested by the current calculation.
    integer(i32), intent(in) :: band_index
    !> First band stored in the evGW0 checkpoint.
    integer(i32), intent(in) :: first_qp_band
    !> Last band stored in the evGW0 checkpoint.
    integer(i32), intent(in) :: last_qp_band

    flag = band_index >= first_qp_band .and. band_index <= last_qp_band
  end function has_checkpoint_qp_energy_for_band

  !> Return whether `eqpsolver` is supported for the initial full explicit evGW0 iteration.
  pure logical function evgw0_seed_eqpsolver_is_supported(eqpsolver)
    !> Quasiparticle-equation solver selected in the GW input.
    integer(i32), intent(in) :: eqpsolver

    select case (eqpsolver)
      case (0, 1, 2)
        evgw0_seed_eqpsolver_is_supported = .true.
      case default
        evgw0_seed_eqpsolver_is_supported = .false.
    end select
  end function evgw0_seed_eqpsolver_is_supported

  !> Return whether a reused explicit evGW0 iteration should warn about ignoring `eqpsolver`.
  pure logical function should_warn_about_ignored_evgw0_eqpsolver(eqpsolver, use_input_qp)
    !> Quasiparticle-equation solver selected in the GW input.
    integer(i32), intent(in) :: eqpsolver
    !> True when the current explicit evGW0 iteration reuses checkpoint quasiparticle energies.
    logical, intent(in) :: use_input_qp

    should_warn_about_ignored_evgw0_eqpsolver = use_input_qp .and. eqpsolver /= 0
  end function should_warn_about_ignored_evgw0_eqpsolver

  !> Assert that the restart files required by the selected split `taskGroup` evGW0 job are present.
  subroutine assert_task_group_restart_files_exist(output_format, is_task_sigmac, is_task_qpeigenvalues, qpoint_vectors, &
      qpeigenvalues_kpoint_indexes, context)
    !> Task-group output format.
    character(len=*), intent(in) :: output_format
    !> True when the current task group includes `sigmac`.
    logical, intent(in) :: is_task_sigmac
    !> True when the current task group includes `QPEigenvalues`.
    logical, intent(in) :: is_task_qpeigenvalues
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> Irreducible k-point indexes used by `QPEigenvalues`.
    integer(i32), intent(in) :: qpeigenvalues_kpoint_indexes(:)
    !> Context string shown in the error message.
    character(len=*), intent(in) :: context

    character(len=40) :: missing_file
    logical :: missing_found

    call get_missing_task_group_restart_file(output_format=output_format, is_task_sigmac=is_task_sigmac, &
      is_task_qpeigenvalues=is_task_qpeigenvalues, qpoint_vectors=qpoint_vectors, &
      qpeigenvalues_kpoint_indexes=qpeigenvalues_kpoint_indexes, missing_file=missing_file, &
      missing_found=missing_found)
    call terminate_if_false(.not. missing_found, trim(context) // ' requires the existing file ' // &
      trim(missing_file) // '.')
  end subroutine assert_task_group_restart_files_exist

  !> Return whether the restart files required by the selected split `taskGroup` evGW0 job are present.
  logical function task_group_restart_files_exist(output_format, is_task_sigmac, is_task_qpeigenvalues, qpoint_vectors, &
      qpeigenvalues_kpoint_indexes, directory) result(flag)
    !> Task-group output format.
    character(len=*), intent(in) :: output_format
    !> True when the current task group includes `sigmac`.
    logical, intent(in) :: is_task_sigmac
    !> True when the current task group includes `QPEigenvalues`.
    logical, intent(in) :: is_task_qpeigenvalues
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> Irreducible k-point indexes used by `QPEigenvalues`.
    integer(i32), intent(in) :: qpeigenvalues_kpoint_indexes(:)
    !> Optional directory that contains the restart files.
    character(len=*), intent(in), optional :: directory

    character(len=40) :: missing_file
    logical :: missing_found

    call get_missing_task_group_restart_file(output_format=output_format, is_task_sigmac=is_task_sigmac, &
      is_task_qpeigenvalues=is_task_qpeigenvalues, qpoint_vectors=qpoint_vectors, &
      qpeigenvalues_kpoint_indexes=qpeigenvalues_kpoint_indexes, missing_file=missing_file, &
      missing_found=missing_found, directory=directory)
    flag = .not. missing_found
  end function task_group_restart_files_exist

  !> Assert that an evGW0 checkpoint matches the current k-point set and QP-band window.
  subroutine assert_checkpoint_matches_expected_qp_window(file_name, kpoint_vectors, first_qp_band, last_qp_band, context)
    !> Checkpoint file to inspect.
    character(len=*), intent(in) :: file_name
    !> Expected irreducible k-point vectors in lattice coordinates.
    real(dp), intent(in) :: kpoint_vectors(:, :)
    !> Expected first QP band.
    integer(i32), intent(in) :: first_qp_band
    !> Expected last QP band.
    integer(i32), intent(in) :: last_qp_band
    !> Context string shown in the error message.
    character(len=*), intent(in) :: context

    call terminate_if_false(checkpoint_matches_expected_qp_window(file_name, kpoint_vectors, first_qp_band, last_qp_band), &
      trim(context) // ' requires checkpoint ' // trim(file_name) // &
      ' to match the current irreducible k-point vectors and QP band window.')
  end subroutine assert_checkpoint_matches_expected_qp_window

  !> Count the restart files required by the selected split `taskGroup` evGW0 job.
  integer(i32) function count_task_group_restart_files(is_task_sigmac, is_task_qpeigenvalues, qpoint_vectors, &
      qpeigenvalues_kpoint_indexes) result(n_files)
    !> True when the current task group includes `sigmac`.
    logical, intent(in) :: is_task_sigmac
    !> True when the current task group includes `QPEigenvalues`.
    logical, intent(in) :: is_task_qpeigenvalues
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> Irreducible k-point indexes used by `QPEigenvalues`.
    integer(i32), intent(in) :: qpeigenvalues_kpoint_indexes(:)

    integer(i32) :: iq

    CALL_ASSERT(size(qpoint_vectors, 1) == 3, 'qpoint_vectors must have size 3 along the first dimension')

    n_files = 0
    if (is_task_sigmac) then
      do iq = 1, size(qpoint_vectors, 2)
        n_files = n_files + 3
        if (gammapoint(qpoint_vectors(:, iq))) n_files = n_files + 3
      end do
    end if
    if (is_task_qpeigenvalues) then
      n_files = n_files + 1 + size(qpeigenvalues_kpoint_indexes)
      if (.not. is_task_sigmac) n_files = n_files + size(qpeigenvalues_kpoint_indexes)
    end if
  end function count_task_group_restart_files

  !> Build the ordered list of restart files required by the selected split `taskGroup` evGW0 job.
  subroutine build_task_group_restart_files(output_format, is_task_sigmac, is_task_qpeigenvalues, qpoint_vectors, &
      qpeigenvalues_kpoint_indexes, file_names)
    !> Task-group output format.
    character(len=*), intent(in) :: output_format
    !> True when the current task group includes `sigmac`.
    logical, intent(in) :: is_task_sigmac
    !> True when the current task group includes `QPEigenvalues`.
    logical, intent(in) :: is_task_qpeigenvalues
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> Irreducible k-point indexes used by `QPEigenvalues`.
    integer(i32), intent(in) :: qpeigenvalues_kpoint_indexes(:)
    !> Restart file names required by the selected tasks.
    character(len=*), intent(out) :: file_names(:)

    character(len=40) :: file_name
    integer(i32) :: ifile, ik, iq
    logical :: file_names_have_expected_size

    file_names_have_expected_size = size(file_names) == count_task_group_restart_files( &
      is_task_sigmac=is_task_sigmac, is_task_qpeigenvalues=is_task_qpeigenvalues, &
      qpoint_vectors=qpoint_vectors, qpeigenvalues_kpoint_indexes=qpeigenvalues_kpoint_indexes)
    CALL_ASSERT(file_names_have_expected_size, 'file_names has the wrong size for the explicit evGW0 task-group restart files')

    ifile = 0
    if (is_task_sigmac) then
      do iq = 1, size(qpoint_vectors, 2)
        ifile = ifile + 1
        call build_file_name('SGI_Q', iq, file_name)
        file_names(ifile) = trim(file_name)

        ifile = ifile + 1
        call build_file_name('BARC_Q', iq, file_name)
        file_names(ifile) = trim(file_name)

        ifile = ifile + 1
        call build_file_name('INVERSE-EPSILON_Q', iq, file_name)
        file_names(ifile) = trim(file_name)

        if (gammapoint(qpoint_vectors(:, iq))) then
          ifile = ifile + 1
          call build_file_name('INVERSE-EPSH', file_name)
          file_names(ifile) = trim(file_name)

          ifile = ifile + 1
          call build_file_name('INVERSE-EPSW1', file_name)
          file_names(ifile) = trim(file_name)

          ifile = ifile + 1
          call build_file_name('INVERSE-EPSW2', file_name)
          file_names(ifile) = trim(file_name)
        end if
      end do
    end if

    if (is_task_qpeigenvalues) then
      ifile = ifile + 1
      file_names(ifile) = get_vxcnn_file_name(output_format)
      do ik = 1, size(qpeigenvalues_kpoint_indexes)
        ifile = ifile + 1
        call build_file_name(file_name_sigmax, qpeigenvalues_kpoint_indexes(ik), file_name)
        file_names(ifile) = trim(file_name)
      end do
      if (.not. is_task_sigmac) then
        do ik = 1, size(qpeigenvalues_kpoint_indexes)
          ifile = ifile + 1
          call build_file_name(file_name_sigmac, qpeigenvalues_kpoint_indexes(ik), file_name)
          file_names(ifile) = trim(file_name)
        end do
      end if
    end if
  end subroutine build_task_group_restart_files

  !> Return the first missing restart file required by the selected split `taskGroup` evGW0 job.
  subroutine get_missing_task_group_restart_file(output_format, is_task_sigmac, is_task_qpeigenvalues, qpoint_vectors, &
      qpeigenvalues_kpoint_indexes, missing_file, missing_found, directory)
    !> Task-group output format.
    character(len=*), intent(in) :: output_format
    !> True when the current task group includes `sigmac`.
    logical, intent(in) :: is_task_sigmac
    !> True when the current task group includes `QPEigenvalues`.
    logical, intent(in) :: is_task_qpeigenvalues
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> Irreducible k-point indexes used by `QPEigenvalues`.
    integer(i32), intent(in) :: qpeigenvalues_kpoint_indexes(:)
    !> First required restart file that is not present.
    character(len=*), intent(out) :: missing_file
    !> True when a required restart file is not present.
    logical, intent(out) :: missing_found
    !> Optional directory that contains the restart files.
    character(len=*), intent(in), optional :: directory

    character(len=40), allocatable :: file_names(:)
    character(len=:), allocatable :: file_path
    integer(i32) :: ifile
    logical :: exists

    allocate(file_names(count_task_group_restart_files(is_task_sigmac=is_task_sigmac, &
      is_task_qpeigenvalues=is_task_qpeigenvalues, qpoint_vectors=qpoint_vectors, &
      qpeigenvalues_kpoint_indexes=qpeigenvalues_kpoint_indexes)))
    call build_task_group_restart_files(output_format=output_format, is_task_sigmac=is_task_sigmac, &
      is_task_qpeigenvalues=is_task_qpeigenvalues, qpoint_vectors=qpoint_vectors, &
      qpeigenvalues_kpoint_indexes=qpeigenvalues_kpoint_indexes, file_names=file_names)

    missing_file = ''
    missing_found = .false.
    do ifile = 1, size(file_names)
      file_path = get_restart_file_path(file_names(ifile), directory)
      inquire(file=trim(file_path), exist=exists)
      if (.not. exists .and. .not. missing_found) then
        missing_file = trim(file_names(ifile))
        missing_found = .true.
      end if
    end do
  end subroutine get_missing_task_group_restart_file

  !> Return the path to a restart file in the current or supplied directory.
  function get_restart_file_path(file_name, directory) result(file_path)
    !> Restart file name.
    character(len=*), intent(in) :: file_name
    !> Optional directory that contains the restart files.
    character(len=*), intent(in), optional :: directory
    !> Path to the restart file.
    character(len=:), allocatable :: file_path

    if (present(directory)) then
      file_path = join_paths(trim(directory), trim(file_name))
    else
      file_path = trim(file_name)
    end if
  end function get_restart_file_path

end module evgw0_validation
