!> Module that contains the types and subroutines needed to execute
!> `taskGroup` in `gw`
module task_group
  use modgw, only: kqset, kset, ibgw, nbgw, kiw, ciw, Gset, Gkset, Gqset, Gqbarc, freq
  use modinput, only: input, gw_type, isspinorb
  use modmpi, only: mpiglobal, terminate_if_false, barrier
  use mod_bands, only: evalfv, numin, nstdf, nstse
  use mod_core_states, only: n_core_states => ncg
  use mod_coulomb_potential, only: calculate_singularities_coeff
  use mod_dielectric_function, only: delete_dielectric_function
  use mod_frequency, only: delete_freqgrid
  use mod_gw_degeneracies, only: delete_degeneracy_module
  use exciting_idiel_interface, only: destroy_idiel_handler
  use mod_kpointset, only: delete_k_vectors, delete_kq_vectors, delete_G_vectors, delete_Gk_vectors
  use mod_selfenergy, only: singc1, singc2, delete_selfenergy
  use precision, only: i32, dp
  use scrcoul_low_dim, only: set_singc12
  use task_Coulomb, only: execute_task_Coulomb
  use task_epsilon, only: execute_task_epsilon
  use task_invertEpsilon, only: execute_task_invertEpsilon
  use task_irreducibleMapping, only: execute_task_irreducibleMapping
  use task_sigmac, only: execute_task_sigmac
  use task_sigmax, only: execute_task_sigmax
  use task_polarizability, only: execute_task_polarizability
  use task_vxc, only: execute_task_vxc
  use task_QPEigenvalues, only: execute_task_QPEigenvalues
  use task_cc4sInterface, only: execute_task_cc4sInterface
  use task_optimized_vxc, only: execute_task_optimized_vxc
  use mod_selfconsistent_gw, only: prepare_next_iteration

  implicit none

  private

  character(len=*), parameter :: task_name = "taskGroup"

  public :: execute_task_group
  public :: task_group_parameters
  public :: initialize_task_group
  public :: execute_task_group_iteration
  public :: deallocate_task_group_global_arrays

  !> Interface to the parameters defined in the input file
  type task_group_parameters
    character(len=20) :: output_format
    logical :: calculate_momentum_matrix
    character(len=20) :: Coulomb_cutoff_type
    character(len=20) :: selfenergy_singularity_treatment
    logical :: task_Coulomb
    logical :: usingIrreducibleWedge_in_task_polarizability = .false.
    logical :: task_polarizability
    logical :: usingIrreducibleWedge_in_task_epsilon = .false.
    logical :: task_epsilon
    logical :: task_invertEpsilon
    logical :: usingIrreducibleWedge_in_task_invertEpsilon = .false.
    logical :: task_irreducibleMapping
    logical :: task_sigmac
    logical :: task_sigmax
    logical :: task_vxc
    logical :: task_optimizedVxc
    logical :: task_prepare_next_selfconsistent_iteration
    logical :: task_QPEigenvalues
    logical :: task_cc4s_interface
    logical :: analytical_limit
    logical :: dry_run
  contains
    procedure :: parse_input
  end type

contains
  !> Execute a group of tasks given in the input file inside the
  !> `taskGroup` element (within `gw`)
  subroutine execute_task_group
    type(task_group_parameters) :: input_parameters

    call input_parameters%parse_input(input%gw)
    call initialize_task_group()
    call execute_task_group_iteration(input_parameters, .true.)
  end subroutine execute_task_group


  !> Execute one `taskGroup` iteration.
  subroutine execute_task_group_iteration(input_parameters, run_full_iteration)
    !> Parsed task-group parameters of the current GW input.
    type(task_group_parameters), intent(in) :: input_parameters
    !> True when the current iteration must run the full task-group setup.
    logical, intent(in) :: run_full_iteration

    integer(i32) :: n_qpoints, n_kpoints, first_state, first_empty_state, last_empty_state, m_dim
    integer(i32) :: n_qpoints_invertepsilon, n_qpoints_epsilon, n_qpoints_coulomb_vertex, i

    n_qpoints = kqset%nkpt
    n_kpoints = kset%nkpt

    call calculate_singularities_coeff(input_parameters%Coulomb_cutoff_type, &
      input_parameters%selfenergy_singularity_treatment, n_qpoints, singc2)

    if (run_full_iteration) then
      if (input_parameters%task_vxc) &
        call execute_task_vxc(ibgw, nbgw, kset%vkl(:, 1:kset%nkpt), input_parameters%output_format, mpiglobal)

      ! clean not used anymore global exciting variables
      call clean_gndstate()

    if( input_parameters%task_Coulomb ) &
      call execute_task_Coulomb( n_qpoints, input_parameters%output_format )

    if( input_parameters%task_cc4s_interface ) then
      call barrier( mpiglobal )

        call terminate_if_false( (n_kpoints==1), &
        'Currently only a single k-point is supported when executing cc4sInterface' )
        call terminate_if_false( (n_qpoints==1), &
        'Currently only a single q-point is supported when executing cc4sInterface' )

      call execute_task_cc4sInterface( input_parameters%output_format  )
    end if

    if( input_parameters%task_sigmax ) then
      ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
      call barrier( mpiglobal )
      call execute_task_sigmax( ibgw, nbgw, n_kpoints, kqset%vqc, input_parameters%output_format )
    end if

      if (input_parameters%task_sigmax) then
        ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
        call barrier(mpiglobal)
        call execute_task_sigmax(ibgw, nbgw, n_kpoints, kqset%vqc, input_parameters%output_format)
      end if

      if (input_parameters%task_polarizability) then
        ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
        call barrier(mpiglobal)

        ! The number of points for the polarizability task depend on the use of symmetry
        ! TODO: There must be a better way than using kset%nkpt
        if (input_parameters%usingIrreducibleWedge_in_task_polarizability) then
          n_qpoints_epsilon = kset%nkpt ! Reduced points
        else
          n_qpoints_epsilon = kqset%nkpt
        end if

        call execute_task_polarizability(n_qpoints_epsilon, input_parameters%output_format)
      end if

      first_empty_state = numin
      last_empty_state = nstdf
      if (input_parameters%task_epsilon) then
        ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
        call barrier(mpiglobal)
        if (input_parameters%usingIrreducibleWedge_in_task_epsilon) then
          call execute_task_epsilon(kqset%vqc, kset%ikp2ik(1:kset%nkpt), kqset%nkpt, first_empty_state, &
            last_empty_state, input_parameters%output_format, input_parameters%dry_run)
        else
          call execute_task_epsilon(kqset%vqc, [(i, i = 1, kqset%nkpt)], kqset%nkpt, first_empty_state, &
            last_empty_state, input_parameters%output_format, input_parameters%dry_run)
        end if
      end if

      if (input_parameters%task_invertEpsilon) then
        ! A barrier is necessary to ensure that all processes have completed outputting the dielectric matrix
        call barrier(mpiglobal)

        ! The number of points for the invert epsilon task depend on the use of symmetry
        ! TODO: See previous comment
        if (input_parameters%usingIrreducibleWedge_in_task_invertEpsilon) then
          n_qpoints_invertepsilon = kset%nkpt ! Reduced points
        else
          n_qpoints_invertepsilon = kqset%nkpt
        end if

        call execute_task_invertEpsilon(n_qpoints_invertepsilon, input_parameters%output_format)
      end if

      if (input_parameters%task_irreducibleMapping) then
        ! A barrier is necessary to ensure that all processes have completed outputting the inverse dielectric matrix
        ! in the irreducible wedge
        call barrier(mpiglobal)
        call execute_task_irreducibleMapping(n_kpoints, input_parameters%output_format)
      end if
    end if

    ! Here, the number of empty states may be different from that used to obtain epsilon
    first_state = 1
    last_empty_state = nstse
    m_dim = last_empty_state
    if (input%gw%coreflag == 'all') m_dim = m_dim + n_core_states
    if (input_parameters%task_sigmac) then
      ! A barrier is necessary to ensure that all processes have completed outputting the inverse of the epsilon
      ! in the full BZ
      call barrier(mpiglobal)
      if (input_parameters%analytical_limit) call set_singc12
      call execute_task_sigmac(n_kpoints, kqset%vqc, first_state, first_state + m_dim - 1, input_parameters%output_format)
    end if

    if (input_parameters%task_QPEigenvalues) then
      ! A barrier is necessary to ensure that all processes have completed outputting sigmac
      call barrier(mpiglobal)
      call execute_task_QPEigenvalues(ibgw, nbgw, kset, input_parameters%output_format)
    end if

    if( input_parameters%task_optimizedVxc ) then
      ! A barrier is necessary to ensure that all processes have completed before starting
      ! the computation of the optimized potential
      call barrier( mpiglobal )
      call execute_task_optimized_vxc( n_kpoints, input_parameters%output_format )
    end if

    if ( input_parameters%task_prepare_next_selfconsistent_iteration ) then
      call barrier( mpiglobal )
      call prepare_next_iteration( input_parameters%output_format )
    end if 

    call delete_selfenergy()
  end subroutine execute_task_group_iteration


  !> Initialize global parameters needed for a GW calculation
  subroutine initialize_task_group()
    ! prepare GW global data
    call init_gw()

    call kintw()
    singc1 = 0.0_dp
    singc2 = 0.0_dp
  end subroutine initialize_task_group


  !> Obtain the parameters defined in the input file that are relevant here
  subroutine parse_input(this, gw_inp)
    !> Parsed task-group parameters to update.
    class(task_group_parameters), intent(inout) :: this
    !> type with the variables given in the input file (inside the gw element)
    type(gw_type), intent(in) :: gw_inp

    call terminate_if_false(associated(gw_inp%taskGroup), &
      'Element taskGroup must be present when taskname='//'"'//task_name//'"')
    call terminate_if_false(associated(gw_inp%barecoul), &
      'Element barecoul must be present when taskname='//'"'//task_name//'"')
    call terminate_if_false(.not. isspinorb(), &
      'Spin-polarized calculations are not currently supported with taskname='//'"'//task_name//'"')
    this%output_format = trim(gw_inp%taskGroup%outputFormat)
    this%calculate_momentum_matrix = .not. gw_inp%rpmat !rpmat means "read pmat"
    this%Coulomb_cutoff_type = trim(gw_inp%barecoul%cutofftype)
    this%selfenergy_singularity_treatment = trim(gw_inp%selfenergy%singularity)
    this%task_Coulomb = associated(gw_inp%taskGroup%Coulomb)
    this%task_polarizability = associated(gw_inp%taskGroup%polarizability)
    if (this%task_polarizability) this%usingIrreducibleWedge_in_task_polarizability = &
      gw_inp%taskGroup%polarizability%usingIrreducibleWedge
    this%task_epsilon = associated(gw_inp%taskGroup%epsilon)
    if (this%task_epsilon) this%usingIrreducibleWedge_in_task_epsilon = gw_inp%taskGroup%epsilon%usingIrreducibleWedge
    this%task_invertEpsilon = associated(gw_inp%taskGroup%invertEpsilon)
    if (this%task_invertEpsilon) this%usingIrreducibleWedge_in_task_invertEpsilon = &
      gw_inp%taskGroup%invertEpsilon%usingIrreducibleWedge
    this%task_irreducibleMapping = associated(gw_inp%taskGroup%irreducibleMapping)
    this%task_sigmac = associated(gw_inp%taskGroup%sigmac)
    this%task_sigmax = associated(gw_inp%taskGroup%sigmax)
    this%task_vxc = associated(gw_inp%taskGroup%vxc)
    this%task_QPEigenvalues = associated(gw_inp%taskGroup%QPEigenvalues)
    this%analytical_limit = (trim(gw_inp%scrcoul%averaging) == '2d' .or. &
      trim(gw_inp%scrcoul%averaging) == 'anisotropic-2d')
    this%dry_run = gw_inp%taskGroup%dryRun
    this%task_cc4s_interface = associated(gw_inp%taskGroup%cc4sInterface)
    this%task_optimizedVxc = associated( gw_inp%taskGroup%optimizedVxc )
    this%task_prepare_next_selfconsistent_iteration = associated( gw_inp%taskGroup%prepareNextSelfconsistentIteration)

    call deallocate_task_group_global_arrays()
  end subroutine parse_input


  !> Deallocate global arrays needed by `taskGroup`
  subroutine deallocate_task_group_global_arrays()
    call delete_dielectric_function(Gamma=.true.)
    call destroy_idiel_handler()
    call delete_degeneracy_module()
    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)
    if (allocated(evalfv)) deallocate(evalfv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
  end subroutine deallocate_task_group_global_arrays
end module task_group
