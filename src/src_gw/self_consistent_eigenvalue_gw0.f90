!> Manage state for the explicit self-consistent eigenvalue GW0 workflow.
!>
!> This module centralizes the iteration bookkeeping introduced by the `gw/evGW0`
!> input element:
!>
!> - discovery of existing checkpoint files `EVALQP_EVGW0xx.OUT`
!> - interpretation of the `do` policy `fromscratch`, `fromfile`, or `auto`
!> - storage of quasiparticle energies reused by later iterations
!> - calculation of the maximum quasiparticle-energy change between iterations
!> - control over whether the current job may write the next workflow checkpoint
!>
!> The explicit workflow supports split HPC execution:
!>
!> - `maxIterations` is the total evGW0 target iteration number, not the number of
!>   loop passes performed inside one job
!> - a restarted split job executes only the current evGW0 iteration and then stops
!>   unless it owns a full irreducible-k `QPEigenvalues` update
!> - only jobs with `QPEigenvalues` on the full irreducible k-point set may write
!>   the next `EVALQP_EVGW0xx.OUT` checkpoint and drive the outer evGW0 loop forward
!>
!> This prevents partial task or subset jobs from silently advancing the workflow
!> with incomplete quasiparticle data.
module self_consistent_eigenvalue_gw0
  use evgw0_validation, only: assert_checkpoint_matches_expected_qp_window, evgw0_seed_eqpsolver_is_supported
  use modinput, only: gw_type
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32, str_32
  use quasiparticle_energies, only: checkpoint_matches_expected_qp_window
  implicit none
  private

  !> Constants
  !> Limit for the number of self-consistent evGW0 iterations
  integer(i32), parameter :: MAX_ITERATIONS = 99
  !> Filenames should have same length with leading zeros
  !> according to the number of digits of MAX_ITERATIONS
  integer(i32), parameter :: NUM_DIGITS = ceiling(log10(real(MAX_ITERATIONS + 1, dp)))

  enum, bind(C)
    enumerator :: evgw0_do_mode
    enumerator :: evgw0_fromscratch, evgw0_fromfile, evgw0_auto
  end enum

  !> Interface to the parameters defined in the `gw/evGW0` input element.
  type evgw0_input_parameters
    private
    logical :: active = .false.
    integer(kind(evgw0_do_mode)) :: do_mode = evgw0_fromscratch
    integer(i32) :: maximum_iterations = 0
    real(dp) :: convergence_tolerance = 0.0_dp
    character(len=20) :: do_string = ''
  contains
    procedure :: get_do_string => evgw0_input_get_do_string
    procedure :: get_max_iterations => evgw0_input_get_max_iterations
    procedure :: get_tolerance => evgw0_input_get_tolerance
    procedure :: is_active => evgw0_input_is_active
    procedure :: is_auto => evgw0_input_is_auto
    procedure :: is_fromfile => evgw0_input_is_fromfile
    procedure :: is_fromscratch => evgw0_input_is_fromscratch
    procedure :: parse_input => parse_evgw0_input
  end type evgw0_input_parameters

  !> Quasiparticle eigenvalues from previous exciting run
  real(dp), allocatable, target :: evalqp_evgw0(:,:)
  !> Status of the current exciting run
  logical :: is_evgw0_active = .false.
  !> True when the explicit evGW0 workflow mode is active
  logical :: is_evgw0_workflow_active = .false.
  !> True when the current iteration uses quasiparticle energies from a previous iteration
  logical :: use_evgw0_input_qp_state = .false.
  !> Current number of self-consistent evGW0 iteration
  integer(i32) :: iteration_number = 0
  !> Maximum absolute quasiparticle-energy difference to the previous iteration
  real(dp) :: evgw0_iteration_max_diff = huge(0.0_dp)
  !> Format string for creating EVALQP filenames that include the iteration number
  character(len=str_32) :: format_str
  !> One-shot override used when `do="auto"` must fall back to a fresh evGW0 start.
  logical :: force_fresh_auto_restart_variable = .false.
  !> True when the current evGW0 job may write the workflow checkpoint for this iteration.
  logical :: write_evgw0_output_checkpoint = .true.

  !> Public types and enumerators
  public :: evgw0_auto
  public :: evgw0_do_mode
  public :: evgw0_fromfile
  public :: evgw0_fromscratch
  public :: evgw0_input_parameters

  !> Public procedures
  public :: finalize_evgw0_iteration
  public :: force_fresh_auto_restart
  public :: get_current_evgw0_input_evalqp_filename
  public :: get_current_evgw0_iteration_number
  public :: get_current_evgw0_output_evalqp_filename
  public :: get_evalqp_evgw0_pointer
  public :: get_evgw0_iteration_max_diff
  public :: has_evgw0_checkpoint
  public :: init_evgw0
  public :: is_evgw0
  public :: is_evgw0_workflow
  public :: pack_evalqp_evgw0_for_kpoints
  public :: reset_evgw0_state
  public :: set_evgw0_output_checkpoint_enabled
  public :: should_write_evgw0_output_checkpoint
  public :: use_evgw0_input_qp
  public :: will_use_evgw0_input_qp

contains

  !> Convert an input `do` string to the internal evGW0 mode enum.
  pure function string_to_evgw0_do_mode(string) result(mode)
    !> Input `do` attribute value.
    character(len=*), intent(in) :: string
    integer(kind(evgw0_do_mode)) :: mode

    select case (trim(string))
      case ('auto')
        mode = evgw0_auto
      case ('fromfile')
        mode = evgw0_fromfile
      case ('fromscratch')
        mode = evgw0_fromscratch
      case default
        mode = evgw0_do_mode
    end select
  end function string_to_evgw0_do_mode

  !> Parse the `gw/evGW0` input element into a small internal interface type.
  subroutine parse_evgw0_input(this, gw_inp)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(inout) :: this
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp

    this%active = associated(gw_inp%evGW0)
    if (this%active) then
      this%do_string = trim(gw_inp%evGW0%do)
      this%do_mode = string_to_evgw0_do_mode(this%do_string)
      this%maximum_iterations = gw_inp%evGW0%maxIterations
      this%convergence_tolerance = gw_inp%evGW0%tolerance
      call terminate_if_false(this%do_mode /= evgw0_do_mode, &
        'Unknown evGW0 do="'//trim(this%do_string)//'"')
      call terminate_if_false(this%maximum_iterations > 0 .and. this%maximum_iterations <= MAX_ITERATIONS, &
        'Error: evGW0 maxIterations must be between 1 and 99')
    else
      this%do_string = ''
      this%do_mode = evgw0_fromscratch
      this%maximum_iterations = 0
      this%convergence_tolerance = 0.0_dp
    end if
  end subroutine parse_evgw0_input

  !> Return the parsed input `do` string.
  pure function evgw0_input_get_do_string(this) result(string)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this
    character(len=20) :: string

    string = this%do_string
  end function evgw0_input_get_do_string

  !> Return the parsed maximum target iteration.
  pure function evgw0_input_get_max_iterations(this) result(max_iterations)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this
    integer(i32) :: max_iterations

    max_iterations = this%maximum_iterations
  end function evgw0_input_get_max_iterations

  !> Return the parsed convergence tolerance.
  pure function evgw0_input_get_tolerance(this) result(tolerance)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this
    real(dp) :: tolerance

    tolerance = this%convergence_tolerance
  end function evgw0_input_get_tolerance

  !> Return whether the `gw/evGW0` input element is present.
  pure logical function evgw0_input_is_active(this) result(flag)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this

    flag = this%active
  end function evgw0_input_is_active

  !> Return whether the input requests automatic restart handling.
  logical function evgw0_input_is_auto(this) result(flag)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this

    flag = this%active .and. this%do_mode == evgw0_auto .and. .not. force_fresh_auto_restart_variable
  end function evgw0_input_is_auto

  !> Return whether the input requires resuming an explicit evGW0 run.
  pure logical function evgw0_input_is_fromfile(this) result(flag)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this

    flag = this%active .and. this%do_mode == evgw0_fromfile
  end function evgw0_input_is_fromfile

  !> Return whether the input requests a fresh explicit evGW0 start.
  logical function evgw0_input_is_fromscratch(this) result(flag)
    !> Parsed evGW0 input parameters.
    class(evgw0_input_parameters), intent(in) :: this

    flag = this%active .and. (this%do_mode == evgw0_fromscratch .or. force_fresh_auto_restart_variable)
  end function evgw0_input_is_fromscratch

  !> Return whether the resolved GW band window covers all quasiparticle states required by explicit evGW0.
  logical function uses_all_qp_states(gw_inp, ib, nb)
    use mod_charge_and_moment, only: chgval

    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> First quasiparticle band used in the current GW calculation.
    integer(i32), intent(in) :: ib
    !> Last quasiparticle band used in the current GW calculation.
    integer(i32), intent(in) :: nb

    integer(i32) :: expected_last_band

    expected_last_band = int(chgval / 2.0_dp, kind=i32) + gw_inp%nempty + 1
    uses_all_qp_states = (ib == 1_i32) .and. (nb == expected_last_band)
  end function uses_all_qp_states

  !> Initialize the filename format used for explicit evGW0 checkpoints.
  subroutine initialize_evgw0_filename_format()
    write(format_str, '(A,I0,A,I0,A)')  &
         '("EVALQP_EVGW0", I', NUM_DIGITS, '.', NUM_DIGITS, ', ".OUT")'
  end subroutine initialize_evgw0_filename_format

  !> Return the highest completed explicit evGW0 checkpoint iteration in the current directory.
  function get_latest_evgw0_checkpoint_iteration() result(iteration)
    integer(i32) :: iteration
    integer(i32) :: icount, io_status, un
    character(len=50) :: filename

    call initialize_evgw0_filename_format()

    iteration = 0
    do icount = MAX_ITERATIONS, 1, -1
       write(filename, format_str) icount
       open(newunit=un, file=filename, status='old', action='read', iostat=io_status)
       if (io_status == 0) close(un)
       if (io_status == 0) then
          iteration = icount
          exit
       end if
    end do
  end function get_latest_evgw0_checkpoint_iteration

  !> Return whether an explicit evGW0 checkpoint exists in the current directory.
  function has_evgw0_checkpoint() result(flag)
    !> True when at least one file `EVALQP_EVGW0xx.OUT` is present.
    logical :: flag

    flag = get_latest_evgw0_checkpoint_iteration() > 0
  end function has_evgw0_checkpoint

  !> True when the next iteration will read quasiparticle energies from a checkpoint.
  logical function will_use_evgw0_input_qp(gw_inp) result(flag)
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    type(evgw0_input_parameters) :: evgw0_input
    integer(i32) :: latest_checkpoint_iteration

    call evgw0_input%parse_input(gw_inp)
    if (.not. evgw0_input%is_active()) then
       flag = .false.
    else if (is_evgw0_workflow_active) then
       flag = .true.
    else
       latest_checkpoint_iteration = get_latest_evgw0_checkpoint_iteration()
       if (evgw0_input%is_fromscratch()) then
          flag = .false.
       else if (evgw0_input%is_auto()) then
          flag = latest_checkpoint_iteration > 0
       else if (evgw0_input%is_fromfile()) then
          call terminate_if_false(latest_checkpoint_iteration > 0, &
               'evGW0 do="fromfile" requires an existing EVALQP_EVGW0 checkpoint file')
          flag = .true.
       else
          call terminate_if_false(.false., 'Unknown evGW0 do="'//trim(evgw0_input%get_do_string())//'"')
          flag = .false.
       end if
    end if
  end function will_use_evgw0_input_qp

  !> Reset all module state associated with the explicit evGW0 workflow.
  !>
  !> This releases stored quasiparticle energies and clears the current iteration
  !> number, `do` mode state, and cached filename format.
  subroutine reset_evgw0_state()
    if (allocated(evalqp_evgw0)) deallocate(evalqp_evgw0)
    is_evgw0_active = .false.
    is_evgw0_workflow_active = .false.
    use_evgw0_input_qp_state = .false.
    iteration_number = 0
    evgw0_iteration_max_diff = huge(0.0_dp)
    format_str = ''
    write_evgw0_output_checkpoint = .true.
  end subroutine reset_evgw0_state

  !> Return whether the current iteration may write an evGW0 workflow checkpoint.
  logical function should_write_evgw0_output_checkpoint() result(flag)
    flag = write_evgw0_output_checkpoint
  end function should_write_evgw0_output_checkpoint

  !> Enable or disable writing the evGW0 workflow checkpoint for the current iteration.
  subroutine set_evgw0_output_checkpoint_enabled(flag)
    !> True when the current job may write `EVALQP_EVGW0xx.OUT`.
    logical, intent(in) :: flag

    write_evgw0_output_checkpoint = flag
  end subroutine set_evgw0_output_checkpoint_enabled

  !> Force the next `do="auto"` initialization to start fresh.
  subroutine force_fresh_auto_restart()
    force_fresh_auto_restart_variable = .true.
  end subroutine force_fresh_auto_restart

  !> Initialize explicit evGW0 state and, when required, read checkpoint quasiparticle energies.
  !>
  !> For the explicit `gw/evGW0` workflow this routine selects the starting
  !> iteration from the `do` policy and, on resumed iterations, reads the
  !> quasiparticle energies from the latest available checkpoint file. The
  !> first full explicit evGW0 iteration may seed the quasiparticle energies with
  !> `eqpsolver=0`, `1`, or `2`; later reuse-based iterations always consume the
  !> stored checkpoint quasiparticle energies instead.
  subroutine init_evgw0(gw_inp, evgw0_input, kset, ib, nb)
    use mod_kpointset, only: k_set
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> Parsed `gw/evGW0` input parameters.
    type(evgw0_input_parameters), intent(in) :: evgw0_input
    !> Irreducible k-point set of the current GW calculation.
    type(k_set),  intent(in)  :: kset
    !> First quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in)  :: ib
    !> Last quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in)  :: nb

    if (.not. evgw0_input%is_active()) then
       call deactivate_evgw0()
    else
       call assert_evgw0_context_is_supported(gw_inp, ib, nb)
       if (is_evgw0_workflow_active) then
          call advance_evgw0_workflow_iteration()
       else
          call initialize_first_evgw0_workflow_iteration(gw_inp, evgw0_input, kset, ib, nb)
       end if
    end if
  end subroutine init_evgw0

  !> Deactivate explicit evGW0 state when the input element is absent.
  subroutine deactivate_evgw0()
    call reset_evgw0_state()
    force_fresh_auto_restart_variable = .false.
  end subroutine deactivate_evgw0

  !> Assert that the current GW setup is supported by explicit evGW0.
  subroutine assert_evgw0_context_is_supported(gw_inp, ib, nb)
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> First quasiparticle band used in the current GW calculation.
    integer(i32), intent(in) :: ib
    !> Last quasiparticle band used in the current GW calculation.
    integer(i32), intent(in) :: nb

    call terminate_if_false(uses_all_qp_states(gw_inp, ib, nb), &
         "Error: self-consistent eigenvalue GW0 requires the full quasiparticle-state window")
    call terminate_if_false(gw_inp%selfenergy%eshift == 0, &
         "Error: self-consistent eigenvalue GW0 does not support alignment of the chemical potential")
    call terminate_if_false(gw_inp%taskname == 'taskGroup', &
         "Error: self-consistent eigenvalue GW0 works currently only with input%gw%taskname='taskGroup'")
    call terminate_if_false(trim(gw_inp%selfenergy%method) == 'ac', &
         "Error: self-consistent eigenvalue GW0 currently supports only input%gw%selfenergy%method='ac'")
  end subroutine assert_evgw0_context_is_supported

  !> Initialize the explicit evGW0 workflow state for the first iteration in this job.
  subroutine initialize_first_evgw0_workflow_iteration(gw_inp, evgw0_input, kset, ib, nb)
    use mod_kpointset, only: k_set
    !> GW input parameters.
    type(gw_type), intent(in) :: gw_inp
    !> Parsed `gw/evGW0` input parameters.
    type(evgw0_input_parameters), intent(in) :: evgw0_input
    !> Irreducible k-point set of the current GW calculation.
    type(k_set), intent(in) :: kset
    !> First quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: ib
    !> Last quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: nb

    integer(i32) :: latest_checkpoint_iteration

    call reset_evgw0_state()
    call initialize_evgw0_filename_format()
    is_evgw0_active = .true.
    is_evgw0_workflow_active = .true.

    latest_checkpoint_iteration = get_latest_evgw0_checkpoint_iteration()
    if (evgw0_input%is_fromscratch()) then
       iteration_number = 1
       use_evgw0_input_qp_state = .false.
    else if (evgw0_input%is_auto()) then
       iteration_number = latest_checkpoint_iteration + 1
       use_evgw0_input_qp_state = latest_checkpoint_iteration > 0
    else if (evgw0_input%is_fromfile()) then
       call terminate_if_false(latest_checkpoint_iteration > 0, &
            'evGW0 do="fromfile" requires an existing EVALQP_EVGW0 checkpoint file')
       iteration_number = latest_checkpoint_iteration + 1
       use_evgw0_input_qp_state = .true.
    else
       call terminate_if_false(.false., 'Unknown evGW0 do="'//trim(evgw0_input%get_do_string())//'"')
    end if

    call terminate_if_false(iteration_number <= MAX_ITERATIONS, &
         "Error: Iteration limit reached. Cannot proceed beyond the defined limit.")

    if (.not. use_evgw0_input_qp_state) then
       call terminate_if_false(evgw0_seed_eqpsolver_is_supported(gw_inp%selfenergy%eqpsolver), &
            "Error: explicit evGW0 currently supports only input%gw%selfenergy%eqpsolver=0, 1, or 2")
    end if

    force_fresh_auto_restart_variable = .false.

    if (use_evgw0_input_qp_state) call read_evgw0_input_checkpoint(evgw0_input, kset, ib, nb)
  end subroutine initialize_first_evgw0_workflow_iteration

  !> Read quasiparticle energies from the current explicit evGW0 input checkpoint.
  subroutine read_evgw0_input_checkpoint(evgw0_input, kset, ib, nb)
    use mod_kpointset, only: k_set
    !> Parsed `gw/evGW0` input parameters.
    type(evgw0_input_parameters), intent(in) :: evgw0_input
    !> Irreducible k-point set of the current GW calculation.
    type(k_set), intent(in) :: kset
    !> First quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: ib
    !> Last quasiparticle band stored in the checkpoint file.
    integer(i32), intent(in) :: nb

    real(dp) :: evalks_temp(ib:nb, kset%nkpt), eferks_temp, eferqp_temp
    character(len=50) :: filename
    logical :: checkpoint_is_usable

    write(filename, format_str) iteration_number - 1
    if (evgw0_input%is_auto()) then
       checkpoint_is_usable = checkpoint_matches_expected_qp_window(trim(filename), kset%vkl(:, 1:kset%nkpt), ib, nb)
    else
       call assert_checkpoint_matches_expected_qp_window(trim(filename), kset%vkl(:, 1:kset%nkpt), ib, nb, &
            'explicit evGW0 restart')
       checkpoint_is_usable = .true.
    end if

    if (checkpoint_is_usable) then
       allocate(evalqp_evgw0(ib:nb, kset%nkpt))
       call readevalqp(trim(filename), kset, ib, nb, evalks_temp, eferks_temp, evalqp_evgw0, eferqp_temp)
    end if
  end subroutine read_evgw0_input_checkpoint

  !> Advance an already active explicit evGW0 workflow to the next iteration.
  subroutine advance_evgw0_workflow_iteration()
    call terminate_if_false(is_evgw0_active, &
         "Error: inconsistent explicit evGW0 state before starting a new iteration.")
    iteration_number = iteration_number + 1
    call terminate_if_false(iteration_number <= MAX_ITERATIONS, &
         "Error: Iteration limit reached. Cannot proceed beyond the defined limit.")
    use_evgw0_input_qp_state = .true.
    call terminate_if_false(allocated(evalqp_evgw0), &
         "Error: explicit evGW0 requires quasiparticle energies from the previous iteration.")
  end subroutine advance_evgw0_workflow_iteration

  !> Return a pointer alias to the quasiparticle energies reused by the current evGW0 iteration.
  !>
  !> Callers use this pointer to avoid copying the full quasiparticle-energy array
  !> in hot paths such as the correlation self-energy evaluation.
  subroutine get_evalqp_evgw0_pointer(array)
    !> Pointer alias to the stored checkpoint or previous-iteration quasiparticle energies.
    real(dp), pointer, intent(out) :: array(:, :)

    call terminate_if_false(is_evgw0_active .and. allocated(evalqp_evgw0), &
         "Error: evalqp_evgw0 array has not been initialized or evGW0 is not enabled.")

    array => evalqp_evgw0
  end subroutine get_evalqp_evgw0_pointer

  !> Return whether the explicit self-consistent eigenvalue GW0 mode is currently active.
  function is_evgw0() result(flag)
    !> True when the explicit evGW0 workflow is active.
    logical :: flag
    flag = is_evgw0_active
  end function is_evgw0

  !> Return whether the explicit multi-iteration evGW0 workflow mode is active.
  function is_evgw0_workflow() result(flag)
    !> True when the explicit `evGW0` input element controls the current run.
    logical :: flag
    flag = is_evgw0_workflow_active
  end function is_evgw0_workflow

  !> Return whether the current iteration reads quasiparticle energies from an earlier iteration.
  function use_evgw0_input_qp() result(flag)
    !> True when the current iteration uses a checkpoint or previous-iteration `EVALQP`.
    logical :: flag
    flag = use_evgw0_input_qp_state
  end function use_evgw0_input_qp

  !> Restrict stored quasiparticle energies to the selected irreducible k-points.
  !>
  !> This is used by tasks that only operate on a selected subset of irreducible
  !> k-points while still reusing checkpoint data from a full previous iteration.
  subroutine pack_evalqp_evgw0_for_kpoints(kpoint_indexes)
    !> Indexes of the irreducible k-points kept for the current task.
    integer(i32), intent(in) :: kpoint_indexes(:)
    real(dp), allocatable :: packed_evalqp(:,:)

    call terminate_if_false(is_evgw0_active, "evGW0 is not enabled.")
    call terminate_if_false(allocated(evalqp_evgw0), &
         "Error: evalqp_evgw0 array has not been initialized.")

    allocate(packed_evalqp(lbound(evalqp_evgw0, 1):ubound(evalqp_evgw0, 1), size(kpoint_indexes)), &
      source=evalqp_evgw0(:, kpoint_indexes))
    call move_alloc(packed_evalqp, evalqp_evgw0)
  end subroutine pack_evalqp_evgw0_for_kpoints

  !> Finalize one evGW0 iteration and store the new quasiparticle energies for the next step.
  !>
  !> When previous quasiparticle energies were reused, this routine also records
  !> the maximum absolute difference between the old and new values. The stored
  !> difference is then used by `task_group` to test the `evGW0` convergence
  !> criterion.
  subroutine finalize_evgw0_iteration(current_evalqp)
    !> Quasiparticle energies obtained for the current iteration.
    real(dp), intent(in) :: current_evalqp(:, :)

    call terminate_if_false(is_evgw0_active, "evGW0 is not enabled.")
    if (use_evgw0_input_qp_state) then
       call terminate_if_false(allocated(evalqp_evgw0), &
            "Error: evalqp_evgw0 array has not been initialized.")
       evgw0_iteration_max_diff = maxval(abs(current_evalqp - evalqp_evgw0))
    else
       evgw0_iteration_max_diff = huge(0.0_dp)
    end if

    if (allocated(evalqp_evgw0)) deallocate(evalqp_evgw0)
    allocate(evalqp_evgw0(lbound(current_evalqp, 1):ubound(current_evalqp, 1), &
         lbound(current_evalqp, 2):ubound(current_evalqp, 2)))
    evalqp_evgw0 = current_evalqp
  end subroutine finalize_evgw0_iteration

  !> Return the maximum quasiparticle-energy change with respect to the previous evGW0 iteration.
  function get_evgw0_iteration_max_diff() result(max_diff)
    !> Maximum absolute quasiparticle-energy difference.
    real(dp) :: max_diff
    max_diff = evgw0_iteration_max_diff
  end function get_evgw0_iteration_max_diff

  !> Return the current explicit evGW0 iteration number.
  function get_current_evgw0_iteration_number() result(number)
    !> Current explicit evGW0 iteration number.
    integer(i32) :: number
    call terminate_if_false(is_evgw0_active, "evGW0 is not enabled.")
    number = iteration_number
  end function get_current_evgw0_iteration_number

  !> Return the checkpoint filename used as quasiparticle-energy input for the current iteration.
  function get_current_evgw0_input_evalqp_filename() result(filename)
    !> Input checkpoint filename for the current iteration.
    character(len=50) :: filename

    ! Ensure evGW0 is enabled
    call terminate_if_false(is_evgw0_active, "evGW0 is not enabled.")

    ! Determine the filename based on the iteration number
    write(filename, format_str) iteration_number - 1
  end function get_current_evgw0_input_evalqp_filename

  !> Return the checkpoint filename written by the current iteration.
  function get_current_evgw0_output_evalqp_filename() result(filename)
    !> Output checkpoint filename for the current iteration.
    character(len=50) :: filename

    ! Ensure evGW0 is enabled
    call terminate_if_false(is_evgw0_active, "evGW0 is not enabled.")

    ! Determine the output filename for the current iteration
    write(filename, format_str) iteration_number
  end function get_current_evgw0_output_evalqp_filename

end module self_consistent_eigenvalue_gw0
