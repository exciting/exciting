!> Module for MPI_comm_split wrapper.
module mod_mpi_comm_split
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  use mpi_f08, only: mpi_comm_split, mpi_comm, mpi_comm_size, mpi_comm_rank, MPI_COMM_NULL
#endif
  use precision, only: i32

  implicit none
  private

  public :: xmpi_comm_split

contains

  !> Split the MPI communicator `parent`.   
  !> All ranks that pass the same `color` will be assigned to the same `child` communicator.
  !> The order of the ranks within `child` will be determined by the order of their `key`.
  subroutine xmpi_comm_split( parent, child, color, key )
    !> parent communicator to split
    type(mpiinfo), intent(inout) :: parent
    !> new communicator
    type(mpiinfo), intent(out) :: child
    !> color
    integer(i32), intent(in) :: color
    !> key
    integer(i32), intent(in) :: key
  
#ifdef MPI
    type(mpi_comm) :: child_comm

    call mpi_comm_split( mpi_comm(parent%comm), color, key, child_comm, parent%ierr )
    child%comm = child_comm%MPI_VAL
    if (child%comm /= MPI_COMM_NULL%MPI_VAL) then
      call mpi_comm_size( child_comm, child%procs, child%ierr )
      call mpi_comm_rank( child_comm, child%rank, child%ierr )
    else
      child%procs = 0
      child%rank = -1
    end if
    child%root = 0
    child%is_root = child%rank == child%root
#else
    child = parent
#endif
  end subroutine xmpi_comm_split

end module mod_mpi_comm_split
