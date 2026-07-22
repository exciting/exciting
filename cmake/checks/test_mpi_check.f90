program test_mpi_check
   use mpi_f08
   implicit none

   integer :: rank, size, ierr

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

   if (size < 1) then
       print *, 'MPI setup error: size is less than 1.'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   else
       print *, 'MPI test passed: Hello from process', rank, 'of', size
   end if

   call MPI_Finalize(ierr)
end program test_mpi_check