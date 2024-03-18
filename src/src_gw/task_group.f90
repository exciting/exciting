!> Module that contains the types and subroutines needed to execute
!> `taskGroup` in `gw`
module task_group
  use modgw, only: kqset
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false
  use mod_coulomb_potential, only: calculate_singularities_coeff
  use mod_selfenergy, only: singc1, singc2
  use precision, only: i32, dp
  use task_Coulomb, only: execute_task_Coulomb
  use task_epsilon, only: execute_task_epsilon

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
    logical :: task_epsilon 
  contains
    procedure :: parse_input
  end type

contains
  !> Execute a group of tasks given in the input file inside the
  !> `taskGroup` element (within `gw`)
  subroutine execute_task_group
    type(task_group_parameters) :: input_parameters
    integer(i32) :: n_qpoints

    call input_parameters%parse_input( input%gw )
    call initialize
    n_qpoints = kqset%nkpt

    call calculate_singularities_coeff( input_parameters%Coulomb_cutoff_type, &
      input_parameters%selfenergy_singularity_treatment, kqset%nkpt, singc2 )

    if( input_parameters%task_Coulomb ) &
      call execute_task_Coulomb( n_qpoints, trim(input_parameters%output_format)=='binary' )

    if( input_parameters%task_epsilon ) &
      call execute_task_epsilon( n_qpoints, input_parameters%output_format )

  end subroutine


  !> Initialize global parameters needed for a GW calculation
  subroutine initialize
    ! prepare GW global data
    call init_gw
      
    ! clean not used anymore global exciting variables
    call clean_gndstate
  
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
      'Element epsilon must be present when taskname='//'"'//task_name//'"' )
    call terminate_if_false( associated(gw_inp%barecoul), &
      'Element barecoul must be present when taskname='//'"'//task_name//'"' )
    this%output_format = trim( gw_inp%taskGroup%outputFormat )
    this%calculate_momentum_matrix = .not. gw_inp%rpmat !rpmat means "read pmat"
    this%Coulomb_cutoff_type = trim( gw_inp%barecoul%cutofftype )
    this%selfenergy_singularity_treatment = trim( gw_inp%selfenergy%singularity )
    this%task_Coulomb = associated( gw_inp%taskGroup%Coulomb )
    this%task_epsilon = associated( gw_inp%taskGroup%epsilon )
  end subroutine
end module