!> Drive the explicit evGW0 workflow while reusing the standard task-group executor.
module evgw0_task_group_driver
  use evgw0_validation, only: assert_checkpoint_matches_expected_qp_window, task_group_restart_files_exist, &
    assert_task_group_restart_files_exist, should_warn_about_ignored_evgw0_eqpsolver
  use gw_info, only: write_to_gwinfo
  use modgw, only: kqset, kset, ibgw, nbgw
  use modinput, only: input
  use mod_kqpts, only: has_full_k_point_coverage, kpoints_sets
  use modmpi, only: mpiglobal, terminate_if_false
  use precision, only: i32, dp
  use quasiparticle_energies, only: checkpoint_matches_expected_qp_window
  use self_consistent_eigenvalue_gw0, only: evgw0_input_parameters, init_evgw0, is_evgw0_workflow, use_evgw0_input_qp, &
    get_evgw0_iteration_max_diff, get_current_evgw0_iteration_number, reset_evgw0_state, &
    get_current_evgw0_input_evalqp_filename, has_evgw0_checkpoint, force_fresh_auto_restart, &
    set_evgw0_output_checkpoint_enabled
  use task_group, only: task_group_parameters, initialize_task_group, execute_task_group_iteration, &
    deallocate_task_group_global_arrays
  use to_char_conversion, only: to_char

  implicit none

  external :: warning

  private

  public :: execute_evgw0_task_group

contains

  !> Execute explicit evGW0 using repeated task-group iterations.
  subroutine execute_evgw0_task_group()
    type(evgw0_input_parameters) :: evgw0_input
    type(task_group_parameters) :: input_parameters

    integer(i32) :: iteration
    real(dp) :: max_diff
    logical :: can_advance_automatically, converged, last_requested_iteration_completed, run_full_iteration, &
      stopped_after_split_job
    character(len=*), parameter :: ignored_eqpsolver_message = &
      'Warning(evGW0): explicit evGW0 uses input selfenergy@eqpsolver only for the initial full iteration. '// &
      'Reuse-based iterations ignore eqpsolver and use checkpoint quasiparticle energies instead.'

    call reset_evgw0_state()
    call evgw0_input%parse_input(input%gw)
    call input_parameters%parse_input(input%gw)

    if (mpiglobal%is_root) then
      call write_to_gwinfo('evGW0 do mode: '//trim(evgw0_input%get_do_string()))
      if (evgw0_input%is_fromscratch() .and. has_evgw0_checkpoint()) then
        call write_to_gwinfo('evGW0 do="fromscratch": ignoring existing checkpoint files')
      end if
    end if

    converged = .false.
    stopped_after_split_job = .false.

    call initialize_first_evgw0_iteration(evgw0_input, input_parameters)
    iteration = get_current_evgw0_iteration_number()
    if (iteration > evgw0_input%get_max_iterations()) then
      call deallocate_task_group_global_arrays()
    else
      call write_evgw0_first_iteration_summary(input_parameters, ignored_eqpsolver_message)
    end if

    do while (iteration <= evgw0_input%get_max_iterations())
      run_full_iteration = .not. use_evgw0_input_qp()
      can_advance_automatically = automatic_evgw0_progress_is_safe(input_parameters, input%gw, kset%nkpt)
      call set_evgw0_output_checkpoint_enabled(can_advance_automatically)

      if (use_evgw0_input_qp()) call validate_evgw0_restart_files(input_parameters, input%gw, &
        kset%vkl(:, 1:kset%nkpt), kqset%vqc, ibgw, nbgw)

      if (mpiglobal%is_root) then
        call write_to_gwinfo('evGW0 iteration '//to_char(iteration))
      end if

      call execute_task_group_iteration_and_cleanup(input_parameters, run_full_iteration)

      if (use_evgw0_input_qp() .and. input_parameters%task_QPEigenvalues .and. can_advance_automatically) then
        max_diff = get_evgw0_iteration_max_diff()
        if (mpiglobal%is_root) then
          call write_to_gwinfo('evGW0 max |dE_QP| = '//to_char(max_diff))
        end if
        if (max_diff < evgw0_input%get_tolerance()) then
          converged = .true.
          exit
        end if
      end if

      if (.not. can_advance_automatically) then
        stopped_after_split_job = .true.
        exit
      end if

      last_requested_iteration_completed = iteration >= evgw0_input%get_max_iterations()
      if (last_requested_iteration_completed) exit

      call initialize_evgw0_iteration(evgw0_input, 'initialization')
      iteration = get_current_evgw0_iteration_number()
    end do

    if (.not. converged .and. mpiglobal%is_root) then
      if (stopped_after_split_job) then
        call write_to_gwinfo('evGW0 stopped after the current job without writing a new workflow checkpoint because automatic progression requires QPEigenvalues on the full irreducible k-point set')
      else
        call write_to_gwinfo('evGW0 reached the maximum number of iterations without satisfying the convergence threshold')
      end if
    end if

    call reset_evgw0_state()
  end subroutine execute_evgw0_task_group

  !> Initialize the first evGW0 task-group iteration and resolve `do="auto"` fallback once.
  subroutine initialize_first_evgw0_iteration(evgw0_input, input_parameters)
    !> Parsed evGW0 input parameters.
    type(evgw0_input_parameters), intent(in) :: evgw0_input
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters

    call initialize_evgw0_iteration(evgw0_input, 'initialization')
    if (evgw0_input%is_auto() .and. use_evgw0_input_qp()) then
      if (.not. auto_restart_evgw0_is_usable(input_parameters, input%gw, kset%vkl(:, 1:kset%nkpt), &
          kqset%vqc, ibgw, nbgw)) then
        if (mpiglobal%is_root) then
          call write_to_gwinfo('evGW0 do="auto": restart prerequisites are incomplete or incompatible')
          call write_to_gwinfo('evGW0 do="auto": falling back to a fresh explicit evGW0 iteration')
        end if
        call reinitialize_auto_restart_from_scratch(evgw0_input)
      end if
    end if
  end subroutine initialize_first_evgw0_iteration

  !> Reinitialize an unusable `do="auto"` restart as a fresh explicit evGW0 iteration.
  subroutine reinitialize_auto_restart_from_scratch(evgw0_input)
    !> Parsed evGW0 input parameters.
    type(evgw0_input_parameters), intent(in) :: evgw0_input

    call reset_evgw0_state()
    call deallocate_task_group_global_arrays()
    call force_fresh_auto_restart()
    call initialize_evgw0_iteration(evgw0_input, 'auto-restart fallback')
  end subroutine reinitialize_auto_restart_from_scratch

  !> Initialize task-group state for the current evGW0 iteration.
  subroutine initialize_evgw0_iteration(evgw0_input, context)
    !> Parsed evGW0 input parameters.
    type(evgw0_input_parameters), intent(in) :: evgw0_input
    !> Context string used in the internal consistency error message.
    character(len=*), intent(in) :: context

    call initialize_task_group()
    call init_evgw0(input%gw, evgw0_input, kset, ibgw, nbgw)
    call terminate_if_false(is_evgw0_workflow(), &
      'Internal error: explicit evGW0 workflow expected after '//trim(context))
  end subroutine initialize_evgw0_iteration

  !> Execute one task-group iteration and release its global arrays before control-flow exits are checked.
  subroutine execute_task_group_iteration_and_cleanup(input_parameters, run_full_iteration)
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters
    !> True when the full task group should run instead of only the selected reused-iteration tasks.
    logical, intent(in) :: run_full_iteration

    call execute_task_group_iteration(input_parameters, run_full_iteration)
    call deallocate_task_group_global_arrays()
  end subroutine execute_task_group_iteration_and_cleanup

  !> Write startup information that applies only to the first evGW0 iteration in this job.
  subroutine write_evgw0_first_iteration_summary(input_parameters, ignored_eqpsolver_message)
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters
    !> Warning issued when `eqpsolver` is ignored by a restarted evGW0 iteration.
    character(len=*), intent(in) :: ignored_eqpsolver_message

    if (mpiglobal%is_root) then
      if (use_evgw0_input_qp()) then
        call write_to_gwinfo('evGW0 restart from '//trim(get_current_evgw0_input_evalqp_filename()))
        call write_evgw0_restart_inputs_summary(input_parameters)
      else
        call write_to_gwinfo('evGW0 starts from iteration 1 without using a checkpoint file')
      end if
      if (should_warn_about_ignored_evgw0_eqpsolver(input%gw%selfenergy%eqpsolver, use_evgw0_input_qp())) then
        call warning(ignored_eqpsolver_message)
        call write_to_gwinfo(ignored_eqpsolver_message)
      end if
    end if
  end subroutine write_evgw0_first_iteration_summary

  !> Return whether the current evGW0 restart inputs are complete enough to resume safely.
  logical function auto_restart_evgw0_is_usable(input_parameters, gw_inp, kpoint_vectors, qpoint_vectors, &
      first_qp_band, last_qp_band) result(flag)
    use modinput, only: gw_type
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> Expected irreducible k-point vectors in lattice coordinates.
    real(dp), intent(in) :: kpoint_vectors(:, :)
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> First quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: first_qp_band
    !> Last quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: last_qp_band

    integer(i32), allocatable :: qp_kpoint_indexes(:)

    flag = checkpoint_matches_expected_qp_window(get_current_evgw0_input_evalqp_filename(), &
      kpoint_vectors, first_qp_band, last_qp_band)
    if (flag) then
      call get_qpeigenvalues_kpoint_indexes(input_parameters, gw_inp, size(kpoint_vectors, 2), qp_kpoint_indexes)
      flag = task_group_restart_files_exist(output_format=input_parameters%output_format, &
        is_task_sigmac=input_parameters%task_sigmac, &
        is_task_qpeigenvalues=input_parameters%task_QPEigenvalues, qpoint_vectors=qpoint_vectors, &
        qpeigenvalues_kpoint_indexes=qp_kpoint_indexes)
    end if
  end function auto_restart_evgw0_is_usable

  !> Validate the restart files required by a reused explicit evGW0 iteration.
  subroutine validate_evgw0_restart_files(input_parameters, gw_inp, kpoint_vectors, qpoint_vectors, &
      first_qp_band, last_qp_band)
    use modinput, only: gw_type
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> Expected irreducible k-point vectors in lattice coordinates.
    real(dp), intent(in) :: kpoint_vectors(:, :)
    !> q-point vectors used to identify restart files.
    real(dp), intent(in) :: qpoint_vectors(:, :)
    !> First quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: first_qp_band
    !> Last quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: last_qp_band

    integer(i32), allocatable :: qp_kpoint_indexes(:)

    call assert_checkpoint_matches_expected_qp_window(get_current_evgw0_input_evalqp_filename(), &
      kpoint_vectors, first_qp_band, last_qp_band, 'explicit evGW0 restart')
    call get_qpeigenvalues_kpoint_indexes(input_parameters, gw_inp, size(kpoint_vectors, 2), qp_kpoint_indexes)
    call assert_task_group_restart_files_exist(output_format=input_parameters%output_format, &
      is_task_sigmac=input_parameters%task_sigmac, &
      is_task_qpeigenvalues=input_parameters%task_QPEigenvalues, qpoint_vectors=qpoint_vectors, &
      qpeigenvalues_kpoint_indexes=qp_kpoint_indexes, context='explicit evGW0 restart')
  end subroutine validate_evgw0_restart_files

  !> Return whether the current job can safely advance the explicit evGW0 loop on its own.
  logical function automatic_evgw0_progress_is_safe(input_parameters, gw_inp, n_kpoints) result(flag)
    use modinput, only: gw_type
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> Total number of irreducible k-points.
    integer(i32), intent(in) :: n_kpoints

    integer(i32), allocatable :: qp_kpoint_indexes(:)

    flag = .false.
    if (input_parameters%task_QPEigenvalues) then
      call get_qpeigenvalues_kpoint_indexes(input_parameters, gw_inp, n_kpoints, qp_kpoint_indexes)
      flag = has_full_k_point_coverage(qp_kpoint_indexes, n_kpoints)
    end if
  end function automatic_evgw0_progress_is_safe

  !> Collect the irreducible k-point indexes used by `QPEigenvalues`.
  subroutine get_qpeigenvalues_kpoint_indexes(input_parameters, gw_inp, n_kpoints, qp_kpoint_indexes)
    use modinput, only: gw_type
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> Total number of irreducible k-points.
    integer(i32), intent(in) :: n_kpoints
    !> Irreducible k-point indexes used by `QPEigenvalues`.
    integer(i32), allocatable, intent(out) :: qp_kpoint_indexes(:)

    type(kpoints_sets) :: qp_k_points

    if (input_parameters%task_QPEigenvalues) then
      call qp_k_points%parse_input(gw_inp%taskGroup%QPEigenvalues%kpointsarray, n_kpoints)
      call qp_k_points%obtain_list_of_indexes()
      allocate(qp_kpoint_indexes(size(qp_k_points%list_of_indexes)))
      qp_kpoint_indexes = qp_k_points%list_of_indexes
    else
      allocate(qp_kpoint_indexes(0))
    end if
  end subroutine get_qpeigenvalues_kpoint_indexes

  !> Write a short summary of the restart files reused by the current evGW0 job.
  subroutine write_evgw0_restart_inputs_summary(input_parameters)
    !> Parsed task-group parameters.
    type(task_group_parameters), intent(in) :: input_parameters

    character(len=160) :: message

    message = 'evGW0 restart reuses existing '
    if (input_parameters%task_sigmac) then
      message = trim(message) // ' SGI, BARC, and inverse-epsilon files'
    end if
    if (input_parameters%task_QPEigenvalues) then
      if (input_parameters%task_sigmac) message = trim(message) // ', plus '
      message = trim(message) // 'VXCNN, SIGMAX, and SIGMAC files'
    end if
    if (.not. input_parameters%task_sigmac .and. .not. input_parameters%task_QPEigenvalues) then
      message = trim(message) // 'checkpoint quasiparticle energies'
    end if
    call write_to_gwinfo(trim(message))
  end subroutine write_evgw0_restart_inputs_summary

end module evgw0_task_group_driver
