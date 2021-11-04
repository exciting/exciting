!> This module only exposes exciting's MPI wrappers,
!> defined in the routines subdirectory, and their
!> corresponding overloads in the serial directory. 
module exciting_mpi
    use mod_mpi_env, only: mpiinfo

#ifdef MPI 
    use mod_mpi_bcast, only: xmpi_bcast
#else
    use mod_serial_bcast, only: xmpi_bcast 
#endif

    implicit none
    public 
    
end module
