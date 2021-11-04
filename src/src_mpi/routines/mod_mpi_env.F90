!> MPI ennvironment type and init. 
!> 
!> Note, any routine in the same module as the MPI environment
!> type (mpiinfo) cannot use exciting's MPI routine bindings, as they 
!> will expect an object of type(mpiinfo) i.e. it creates a circular 
!> dependency.
module mod_mpi_env
#ifdef MPI 
   use mpi
#endif 

   use precision, only: sp

   implicit none
   private

   !> Exposed MPI environment type, and MPI initialisation
   public :: mpiinfo

   !> MPI environment type, containing a description of the MPI environment
   type mpiinfo
      integer(sp) :: comm   !! MPI communicator
      integer(sp) :: rank   !! Rank or process ID
      integer(sp) :: procs  !! Number of MPI processes
      integer(sp) :: ierr   !! Error integer
      integer(sp) :: root   !! Root process associated with comm
      logical :: is_root    !! Indicates if the process is root
   contains
      procedure :: init => xmpi_init
   end type mpiinfo

contains

   !> Initialise the MPI_COMM_WORLD, and the MPI environment object
   !>
   !> TODO(Alex) Issue 26. Extend this to threading support: mpi_init_thread()
   subroutine xmpi_init(mpi_env)
    !> MPI environment 
    class(mpiinfo), intent(inout) :: mpi_env

#ifdef MPI
    call mpi_init(mpi_env%ierr)
    mpi_env%comm = mpi_comm_world
    call mpi_comm_size(mpi_env%comm, mpi_env%procs, mpi_env%ierr)
    call mpi_comm_rank(mpi_env%comm, mpi_env%rank, mpi_env%ierr)
#else
    mpi_env%comm = 0
    mpi_env%procs = 1
    mpi_env%rank= 0
#endif

    ! Root process is always defined as 0
    ! Assumes local indexing of root ids for sub-communicators
    mpi_env%root = 0
    mpi_env%is_root = mpi_env%rank == mpi_env%root

   end subroutine xmpi_init

end module mod_mpi_env
