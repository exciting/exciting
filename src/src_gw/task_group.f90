!> Module that contains the types and subroutines needed to execute
!> `taskGroup` in `gw`
module task_group
  use modgw, only: kqset, kset, ibgw, nbgw, kiw, ciw, Gset, Gkset, Gqset, Gqbarc, freq, nvelgw, nbandsgw
  use modinput, only: input, gw_type, isspinorb
  use modmpi, only: mpiglobal, terminate_if_false, barrier
  use mod_bands, only: evalfv, numin, nstdf, nstse
  use mod_core_states, only: n_core_states => ncg
  use mod_coulomb_potential, only: calculate_singularities_coeff
  use mod_dielectric_function, only: delete_dielectric_function
  use mod_frequency, only: delete_freqgrid
  use mod_kpointset, only: delete_k_vectors, delete_kq_vectors, delete_G_vectors, &
    & delete_Gk_vectors
  use mod_selfenergy, only: singc1, singc2, evalqp, delete_selfenergy
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

  implicit none
  
  private

  character(len=*), parameter :: task_name = "taskGroup"

  public :: execute_task_group

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
    logical :: task_QPEigenvalues
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
    integer(i32) :: n_qpoints, n_kpoints, first_state, first_empty_state, last_empty_state, m_dim
    integer(i32) :: n_qpoints_invertepsilon, n_qpoints_epsilon, i

    call input_parameters%parse_input( input%gw )
    call initialize()
    n_qpoints         = kqset%nkpt
    n_kpoints         = kset%nkpt

    call calculate_singularities_coeff( input_parameters%Coulomb_cutoff_type, &
      input_parameters%selfenergy_singularity_treatment, n_qpoints, singc2 )

    if( input_parameters%task_vxc ) &
      call execute_task_vxc( ibgw, nbgw, kset%vkl(:, 1:kset%nkpt), input_parameters%output_format, mpiglobal )

    ! clean not used anymore global exciting variables
    call clean_gndstate()

    if( input_parameters%task_Coulomb ) &
      call execute_task_Coulomb( n_qpoints, input_parameters%output_format )
    
    if( input_parameters%task_sigmax ) then
      ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
      call barrier( mpiglobal )
      call execute_task_sigmax( ibgw, nbgw, n_kpoints, kqset%vqc, input_parameters%output_format )
    end if

    if (input_parameters%task_polarizability) then
      ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
      call barrier( mpiglobal )
      
      ! The number of points for the polarizability task depend on the use of symmetry
      ! TODO: There must be a better way than using kset%nkpt
      if (input_parameters%usingIrreducibleWedge_in_task_polarizability) then
        n_qpoints_epsilon = kset%nkpt ! Reduced points
      else
        n_qpoints_epsilon = kqset%nkpt
      end if

      call execute_task_polarizability( n_qpoints_epsilon, input_parameters%output_format )
    end if

    first_empty_state = numin
    last_empty_state = nstdf
    if( input_parameters%task_epsilon ) then
      ! A barrier is necessary to ensure that all processes have completed outputting the bare Coulomb matrix
      call barrier( mpiglobal )
      if (input_parameters%usingIrreducibleWedge_in_task_epsilon) then
        call execute_task_epsilon( kqset%vqc, kset%ikp2ik(1:kset%nkpt), kqset%nkpt, first_empty_state, last_empty_state, input_parameters%output_format, input_parameters%dry_run )
      else
        call execute_task_epsilon( kqset%vqc, [(i, i = 1, kqset%nkpt)], kqset%nkpt, first_empty_state, last_empty_state, input_parameters%output_format, input_parameters%dry_run )
      end if
    end if

    if( input_parameters%task_invertEpsilon ) then
      ! A barrier is necessary to ensure that all processes have completed outputting the dielectric matrix
      call barrier( mpiglobal )
      
      ! The number of points for the invert epsilon task depend on the use of symmetry
      ! TODO: See previous comment
      if (input_parameters%usingIrreducibleWedge_in_task_invertEpsilon) then
        n_qpoints_invertepsilon = kset%nkpt ! Reduced points
      else
        n_qpoints_invertepsilon = kqset%nkpt
      end if

      call execute_task_invertEpsilon( n_qpoints_invertepsilon, input_parameters%output_format )
    end if

    if ( input_parameters%task_irreducibleMapping ) then
      ! A barrier is necessary to ensure that all processes have completed outputting the inverse dielectric matrix
      ! in the irreducible wedge
      call barrier( mpiglobal )
      call execute_task_irreducibleMapping(n_kpoints, input_parameters%output_format)
    end if

    ! Here, the number of empty states may be different from that used to obtain epsilon
    first_state = 1
    last_empty_state = nstse
    m_dim = last_empty_state
    if( input%gw%coreflag == 'all' ) m_dim = m_dim + n_core_states
    if( input_parameters%task_sigmac ) then
      ! A barrier is necessary to ensure that all processes have completed outputting the inverse of the epsilon
      ! in the full BZ
      call barrier( mpiglobal )
      if( input_parameters%analytical_limit ) call set_singc12
      call execute_task_sigmac( n_kpoints, kqset%vqc, first_state, first_state+m_dim-1, input_parameters%output_format )
    end if
    
    if( input_parameters%task_QPEigenvalues ) then
      ! A barrier is necessary to ensure that all processes have completed outputting sigmac
      call barrier( mpiglobal )
      call execute_task_QPEigenvalues( ibgw, nbgw, kset, input_parameters%output_format )
    end if

    call delete_selfenergy()

  end subroutine


  !> Initialize global parameters needed for a GW calculation
  subroutine initialize()
    ! prepare GW global data
    call init_gw()
      
    call kintw()
    singc1 = 0.0_dp
    singc2 = 0.0_dp
  end subroutine


  !> Obtain the parameters defined in the input file that are relevant here
  subroutine parse_input( this, gw_inp )
    class(task_group_parameters), intent(inout) :: this
    !> type with the variables given in the input file (inside the gw element)
    type(gw_type), intent(in) :: gw_inp

    call terminate_if_false( associated(gw_inp%taskGroup), &
      'Element taskGroup must be present when taskname='//'"'//task_name//'"' )
    call terminate_if_false( associated(gw_inp%barecoul), &
      'Element barecoul must be present when taskname='//'"'//task_name//'"' )
    call terminate_if_false( .not. isspinorb(), &
      'Spin-polarized calculations are not currently supported with taskname='//'"'//task_name//'"' )
    this%output_format = trim( gw_inp%taskGroup%outputFormat )
    this%calculate_momentum_matrix = .not. gw_inp%rpmat !rpmat means "read pmat"
    this%Coulomb_cutoff_type = trim( gw_inp%barecoul%cutofftype )
    this%selfenergy_singularity_treatment = trim( gw_inp%selfenergy%singularity )
    this%task_Coulomb = associated( gw_inp%taskGroup%Coulomb )
    this%task_polarizability = associated( gw_inp%taskGroup%polarizability)
    if (this%task_polarizability) this%usingIrreducibleWedge_in_task_polarizability = gw_inp%taskGroup%polarizability%usingIrreducibleWedge 
    this%task_epsilon = associated( gw_inp%taskGroup%epsilon )
    if (this%task_epsilon) this%usingIrreducibleWedge_in_task_epsilon = gw_inp%taskGroup%epsilon%usingIrreducibleWedge
    this%task_invertEpsilon = associated( gw_inp%taskGroup%invertEpsilon )
    if (this%task_invertEpsilon) this%usingIrreducibleWedge_in_task_invertEpsilon = gw_inp%taskGroup%invertEpsilon%usingIrreducibleWedge
    this%task_irreducibleMapping = associated( gw_inp%taskGroup%irreducibleMapping )
    this%task_sigmac = associated( gw_inp%taskGroup%sigmac )
    this%task_sigmax = associated( gw_inp%taskGroup%sigmax )
    this%task_vxc = associated( gw_inp%taskGroup%vxc )
    this%task_QPEigenvalues = associated( gw_inp%taskGroup%QPEigenvalues )
    this%analytical_limit = ( trim(gw_inp%scrcoul%averaging) == '2d' .or. trim(gw_inp%scrcoul%averaging) == 'anisotropic-2d')
    this%dry_run = gw_inp%taskGroup%dryRun

    call deallocate_global_arrays

  end subroutine

  !> Deallocate global arrays needed
  subroutine deallocate_global_arrays
    call delete_dielectric_function( Gamma=.true. )
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
  end subroutine
end module
