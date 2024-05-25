!> This module only exposes exciting's MPI wrappers,
!> defined in the routines subdirectory, and their
!> corresponding overloads in the serial directory. 
module exciting_mpi
    use mod_mpi_env, only: mpiinfo
    
#ifdef MPI 
    use mod_mpi_bcast, only: xmpi_bcast
    use mod_mpi_allgather, only: xmpi_allgather, xmpi_allgatherv
    use mod_mpi_allreduce, only: xmpi_allreduce
#else
    use mod_serial_bcast, only: xmpi_bcast
    use mod_serial_allgather, only: xmpi_allgather, xmpi_allgatherv
    use mod_serial_allreduce, only: xmpi_allreduce
#endif

    implicit none
    public 
    
end module
