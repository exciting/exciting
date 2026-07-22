!> This module only exposes exciting's MPI wrappers,
!> defined in the routines subdirectory
module exciting_mpi
  use mod_mpi_env, only: mpiinfo
  use mod_mpi_allgather, only: xmpi_allgather, xmpi_allgatherv
  use mod_mpi_allreduce, only: xmpi_allreduce
  use mod_mpi_bcast, only: xmpi_bcast
  use mod_mpi_gather, only: xmpi_gather, xmpi_gatherv
  use mod_mpi_reduce, only: xmpi_reduce
  use mod_mpi_comm_split, only: xmpi_comm_split
  
  implicit none
  public 
    
end module
