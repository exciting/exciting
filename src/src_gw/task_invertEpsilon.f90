module task_invertEpsilon
  use constants, only: zzero
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use invert_dielectric_function, only: calcinveps
  use modinput, only: input, gw_type
  use modgw, only: freq, kqset, kset, time_dfinv
  use modmpi, only: terminate_if_false, mpiglobal, distribute_loop
  use modxs, only: symt2
  use mod_dielectric_function, only: epsilon, eps00, epsh, epsw1, epsw2, &
    read_epsilon_from_file, write_inverse_epsilon_to_file
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use precision, only: i32, max_length => str_32, str_128

  implicit none

  private

  public :: execute_task_invertEpsilon

  character(len=*), parameter :: task_name = "invertEpsilon"

  !> Interface to the parameters defined in the input file
  type task_invertEpsilon_parameters
    private
    type(kpoints_sets) :: q_points
    integer(i32) :: n_omega
    logical :: usingIrreducibleWedge
    character(len=max_length) :: output_format
  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_qpt )
  class(task_invertEpsilon_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp
  !> maximum number of q-points
  integer(i32),  intent(in) :: n_qpt

  call this%sanity_checks( gw_inp )
  this%usingIrreducibleWedge = gw_inp%taskGroup%invertEpsilon%usingIrreducibleWedge
  call this%q_points%parse_input( gw_inp%taskGroup%invertEpsilon%qpointsarray, n_qpt )
  this%n_omega = gw_inp%freqgrid%nomeg
  this%output_format = trim( gw_inp%taskGroup%outputFormat )
end subroutine

subroutine sanity_checks( this, gw_inp )
  class(task_invertEpsilon_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%invertEpsilon), &
    'Element epsilon must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%freqgrid), &
    'Element freqgrid must be present when executing '//'"'//task_name//'"' )
end subroutine


subroutine execute_task_invertEpsilon( n_qpoints_max, file_format )
  integer(i32), intent(in) :: n_qpoints_max
  character(len=*), intent(in) :: file_format

  integer(i32) :: iq, iq_IO, iq_reducible 
  integer(i32) :: i, i_start, i_end, omega_i, omega_f
  type(task_invertEpsilon_parameters) :: input_parameters
  character(len=str_128) :: string
  type(mpiinfo) :: mpi_environment_qpoints
  
  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  
  call input_parameters%parse_input( input%gw, n_qpoints_max )
  call input_parameters%q_points%obtain_list_of_indexes()
  call mpi_environment_qpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_qpoints, size(input_parameters%q_points%list_of_indexes), &
    i_start, i_end )

  omega_i = 1
  omega_f = input_parameters%n_omega
  do i = i_start, i_end
    iq = input_parameters%q_points%list_of_indexes(i)
    if (input_parameters%usingIrreducibleWedge) then
      iq_IO        = iq ! iq refers to the irreducible wedge
      iq_reducible = kset%ikp2ik(iq)
    else 
      iq_reducible = iq ! iq refers to the full BZ
      iq_IO        = iq_reducible
    end if

    if( mpiglobal%rank == 0 ) then
      if (input_parameters%usingIrreducibleWedge) then
        write( string, * ) '('//task_name//'): q-point cycle, iq (irreducible) = ', iq
      else
        write( string, * ) '('//task_name//'): q-point cycle, iq = ', iq
      end if 
      call write_to_gwinfo( string )
    end if
    
    ! If only working with irreducible points ignore non representative points
    if (input_parameters%usingIrreducibleWedge .and. iq_reducible /= kset%ikp2ik(kset%ik2ikp(iq_reducible))) cycle
    
    Gamma = gammapoint(kqset%vqc(:,iq_reducible))
    call read_epsilon_from_file( iq_IO, Gamma, file_format, input_parameters%usingIrreducibleWedge)
    if( Gamma ) then
      if( .not. allocated(eps00) ) allocate(eps00(3, 3, omega_i:omega_f), source=zzero )
      call calcinveps( omega_i, omega_f, Gamma, input%gw%scrcoul, freq%fconv, symt2,&
                      &epsilon, epsw1, epsw2, epsh, eps00, time_dfinv)
    else
      call calcinveps( omega_i, omega_f, Gamma, freqtype=freq%fconv, epsilon=epsilon, time_dfinv=time_dfinv)
    end if
    call write_inverse_epsilon_to_file( iq_IO, Gamma, file_format, input_parameters%usingIrreducibleWedge)
  end do

end subroutine

end module
