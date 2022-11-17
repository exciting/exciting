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
   use precision, only: i32

   implicit none
   private

   !> Exposed MPI environment type, and MPI initialisation
   public :: mpiinfo

   !> MPI environment type, containing a description of the MPI environment
   type mpiinfo
      integer(i32) :: comm   !! MPI communicator
      integer(i32) :: rank   !! Rank or process ID
      integer(i32) :: procs  !! Number of MPI processes
      integer(i32) :: ierr   !! Error integer
      integer(i32) :: root   !! Root process associated with comm
      logical :: is_root     !! Indicates if the process is root
   contains
      procedure :: init => xmpi_init
   end type mpiinfo

contains

   !> Initialise the MPI_COMM_WORLD, and the MPI environment object
   !>
   !> TODO(Alex) Issue 26. Extend this to threading support: mpi_init_thread()
   subroutine xmpi_init(this, comm)
    !> MPI environment 
    class(mpiinfo), intent(inout) :: this
    !> MPI communicator
    integer, intent(in) :: comm
#ifdef MPI
    this%comm = comm
    call mpi_comm_size(this%comm, this%procs, this%ierr)
    call mpi_comm_rank(this%comm, this%rank, this%ierr)
#else
    this%comm = 0
    this%procs = 1
    this%rank= 0
#endif

    ! Root process is always defined as 0
    ! Assumes local indexing of root ids for sub-communicators
    this%root = 0
    this%is_root = this%rank == this%root

   end subroutine xmpi_init

end module mod_mpi_env
