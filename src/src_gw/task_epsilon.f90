!> This module contains classes and subroutines needed to execute the 
!> task `epsilon`, used as an element of `taskGroup` in `gw`
module task_epsilon
  use asserts, only: assert
  use exciting_mpi, only: mpiinfo
  use gw_io, only: write_to_file, build_file_name, write_to_gwinfo, write_to_gwinfo_boxmessage
  use math_utils, only: all_zero
  use modgw, only: freq, fgw, kset, kqset, Gset, Gqbarc, Gkset, Gqset, ciw, kiw
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, mpiglobal, distribute_loop
  use mod_bands, only: evalfv
  use mod_coulomb_potential, only: read_coulomb_potential_from_file
  use mod_dielectric_function, only: write_epsilon_to_file, init_dielectric_function, delete_dielectric_function
  use mod_frequency, only: delete_freqgrid
  use mod_kpointset, only: delete_k_vectors, delete_kq_vectors, delete_G_vectors, &
    & delete_Gk_vectors
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma
  use mod_product_basis, only: mbsiz, read_sgi_from_file
  use mod_selfenergy, only: singc1, singc2
  use precision, only: dp, i32

  implicit none

  private

  public :: execute_task_epsilon

  integer(i32), parameter :: max_length = 30

  character(len=*), parameter :: task_name = "epsilon"

  !> Interface to the parameters defined in the input file
  type task_epsilon_parameters
    private
    type(kpoints_sets) :: q_points
    integer(i32) :: n_omega
    character(len=max_length) :: output_format
    character(len=max_length) :: screened_coulomb_model
    logical :: calculate_momentum_matrix
  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_qpt )
  class(task_epsilon_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_qpt

  call this%sanity_checks( gw_inp )
  call this%q_points%parse_input( gw_inp%taskGroup%epsilon%qpointsarray, n_qpt )
  this%n_omega = gw_inp%freqgrid%nomeg
  this%screened_coulomb_model = gw_inp%scrcoul%scrtype
  this%output_format = trim( adjustl( gw_inp%taskGroup%outputFormat ) )
  this%calculate_momentum_matrix = .not. gw_inp%rpmat !rpmat means "read pmat"
end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  class(task_epsilon_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%epsilon), &
    'Element epsilon must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%barecoul), &
    'Element barecoul must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%scrcoul), &
    'Element scrcoul must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%freqgrid), &
    'Element freqgrid must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( trim(gw_inp%scrcoul%scrtype)=='rpa', &
    'Only the RPA method is currently supported when executing '//'"'//task_name//'"' )
end subroutine


!> Subroutine to be invoked when task `epsilon` must be executed
subroutine execute_task_epsilon( n_qpoints_max, file_format )
  integer(i32), intent(in) :: n_qpoints_max
  character(len=*), intent(in) :: file_format

  integer(i32) :: iq, i, i_start, i_end, omega_i, omega_f
  integer(i32), parameter :: maxlen=60
  character(len=maxlen) :: string
  type(task_epsilon_parameters) :: input_parameters
  type(mpiinfo) :: mpi_environment_qpoints
  
  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  call input_parameters%parse_input( input%gw, n_qpoints_max )
  call input_parameters%q_points%obtain_list_of_indexes()
  call mpi_environment_qpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_qpoints, size(input_parameters%q_points%list_of_indexes), &
    i_start, i_end )
  
  omega_i = 1
  omega_f = input_parameters%n_omega
  ! Attention: calcpmatgw makes use of MPI parallelization and calls a mpi_barrier
  if( isGammaInList( kqset%vqc(:,input_parameters%q_points%list_of_indexes) ) ) call calcpmatgw
  do i = i_start, i_end
    iq = input_parameters%q_points%list_of_indexes(i)
    if( mpiglobal%rank == 0) then
      write( string, * ) '('//task_name//'): q-point cycle, iq = ', iq
      call write_to_gwinfo( string )
    end if
    call read_coulomb_potential_from_file( iq, file_format )
    call read_sgi_from_file( iq, file_format )
    call calcmpwipw( iq )
    Gamma = all_zero( kqset%vqc(:,iq), tol=1.e-6_dp )
    call init_dielectric_function( mbsiz, omega_i, omega_f, Gamma )
    call calcepsilon( iq, omega_i, omega_f )
    call write_epsilon_to_file( iq, Gamma, file_format )
  end do

  call deallocate_global_arrays
end subroutine


!> Check if the gamma point is among the q-points
pure logical function isGammaInList( q_points )
  !> List of q-points (first dimension has the 3 cartesian coordinates)
  real(dp), intent(in) :: q_points(:, :)

  integer(i32) :: i

  isGammaInList = .false.
  do i = 1, size( q_points, 2 )
    isGammaInList = all_zero( q_points(:, i), tol=1.e-6_dp )
    if( isGammaInList ) exit
  end do
end function


!> Deallocate global arrays needed to obtain the dielectric matrix
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