
subroutine barrier(rank,procs,un,async,string)
  implicit none
  ! arguments
  integer, intent(in) :: rank,procs,un,async
  character(*) :: string
  ! local variables
  character(*), parameter :: thisnam = 'barrier'

  ! do nothing if only one process
  if (procs.eq.1) return

  ! call the MPI barrier
#ifdef MPI
  call MPI_barrier(mpi_comm_world,ierr)
#endif

end subroutine barrier
