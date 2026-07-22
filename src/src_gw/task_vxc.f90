!> Module designed for the task vxc
module task_vxc
#include "asserts.fpp"
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo_boxmessage
  use mod_kqpts, only: kpoints_sets
  use mod_gw_degeneracies, only: ibgw_including_degeneracy, nbgw_including_degeneracy
  use mod_vxc, only: calcvxcnn, write_vxcnn, deallocate_vxcnn
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32

  implicit none 

  private

  character(len=*), parameter :: task_name = "vxc"

  !> Interface to the parameters defined in the input file
  type task_vxc_parameters
    private
    type(kpoints_sets) :: k_points
  contains
    procedure :: parse_input, sanity_checks
  end type

  public :: execute_task_vxc

contains
!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_kpt )
  !> VXC task parameters to update.
  class(task_vxc_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_kpt

  call this%sanity_checks( gw_inp )
  call this%k_points%parse_input( gw_inp%taskGroup%vxc%kpointsarray, n_kpt )

end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  !> VXC task parameters to check.
  class(task_vxc_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%vxc), &
    'Element vxc must be present when executing '//'"'//task_name//'"' )

end subroutine


!> Execute a task vxc calculation
subroutine execute_task_vxc( first_band, last_band, kpt_latt_coord, file_format, mpi_env )
  !> Index of the first KS band for which vxc is evaluated
  integer(i32), intent(in) :: first_band
  !> Index of the last KS band for which vxc is evaluated
  integer(i32), intent(in) :: last_band
  !> Irreducible k-points in lattice coordinates
  real(dp), intent(in) :: kpt_latt_coord(:, :)
  !> Format of input/ouput files
  character(len=*), intent(in) :: file_format
  !> The MPI environment type, for distribution over MPI processes
  type(mpiinfo), intent(in) :: mpi_env

  integer(i32) :: n_kpoints_max
  type(task_vxc_parameters) :: input_parameters
  logical :: my_rank_writes_outputs

  ! kpt_latt_coord should have 3 coordinates for each k-point
  CALL_ASSERT( size( kpt_latt_coord, 1 ) == 3,  'Size of kpt_latt_coord along 1st dim. must be = 3')
  ! first_band, last_band should fit in the range [ibgw_including_degeneracy, nbgw_including_degeneracy]
  CALL_ASSERT( first_band >= ibgw_including_degeneracy, 'first_band must be >= ibgw_including_degeneracy' )
  CALL_ASSERT( last_band <= nbgw_including_degeneracy, 'last_band must be >= nbgw_including_degeneracy' )

  n_kpoints_max = size( kpt_latt_coord, 2 )
  call input_parameters%parse_input( input%gw, n_kpoints_max )
  call input_parameters%k_points%obtain_list_of_indexes()

  my_rank_writes_outputs = ( mpi_env%rank == mpi_env%root )
  if( my_rank_writes_outputs ) call write_to_gwinfo_boxmessage( '=', 'task: ' // task_name )
  associate( list => input_parameters%k_points%list_of_indexes )
    call calcvxcnn( ibgw_including_degeneracy, nbgw_including_degeneracy, list, kpt_latt_coord(:, list), mpi_env )
  end associate
  if( my_rank_writes_outputs ) call write_vxcnn( file_format, first_band, last_band )
  call deallocate_vxcnn
  
end subroutine

end module
