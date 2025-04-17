!> exciting wrappers for mpi_bcast
module mod_mpi_bcast
#ifdef MPI
  use mpi_f08, only: mpi_bcast, mpi_comm, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
    MPI_DOUBLE_COMPLEX, MPI_CHARACTER, MPI_LOGICAL
#endif
  use mod_mpi_env, only: mpiinfo
  use precision, only: dp, i32

  implicit none
  private

  !> exciting wrapper of mpi_bcast 
  public :: xmpi_bcast

  interface xmpi_bcast
      module procedure :: &
        mpi_bcast_logical, &
        mpi_bcast_integer_i32, &
        mpi_bcast_real_dp, &
        mpi_bcast_complex_dp, &
        mpi_bcast_character
  end interface

contains

  !> Broadcasts logical from the process with rank bcasting_rank (root by default) to all other processes of the group.
  subroutine mpi_bcast_logical( mpi_env, buffer, bcasting_rank )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer
    logical, intent(inout) :: buffer(..)
    !> Broadcasting rank
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank
    bc_rank = mpi_env%root
    if ( present ( bcasting_rank ) ) bc_rank = bcasting_rank
    call mpi_bcast( buffer, size( buffer ), MPI_LOGICAL, bc_rank, &
      mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine


  !> Broadcasts integer i32 array from the process with rank bcasting_rank (root by default) to all other processes of the group. 
  subroutine mpi_bcast_integer_i32( mpi_env, buffer, bcasting_rank )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer
    integer(i32), intent(inout) :: buffer(..)
    !> Broadcasting rank
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank
    bc_rank = mpi_env%root
    if ( present ( bcasting_rank ) ) bc_rank = bcasting_rank
    call mpi_bcast( buffer, size( buffer ), MPI_INTEGER, bc_rank, &
      mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine 


  !> Broadcasts real dp array from the process with rank bcasting_rank (root by default) to all other processes of the group. 
  subroutine mpi_bcast_real_dp( mpi_env, buffer, bcasting_rank )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer
    real(dp), intent(inout) :: buffer(..)
    !> Broadcasting rank
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank
    bc_rank = mpi_env%root
    if ( present ( bcasting_rank ) ) bc_rank = bcasting_rank
    call mpi_bcast( buffer, size( buffer ), MPI_DOUBLE_PRECISION, bc_rank, &
      mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine 

  
  !> Broadcasts complex dp array from the process with rank bcasting_rank (root by default) to all other processes of the group. 
  subroutine mpi_bcast_complex_dp( mpi_env, buffer, bcasting_rank )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer
    complex(dp), intent(inout) :: buffer(..)
    !> Broadcasting rank
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank
    bc_rank = mpi_env%root
    if ( present ( bcasting_rank ) ) bc_rank = bcasting_rank
    call mpi_bcast( buffer, size( buffer ), MPI_DOUBLE_COMPLEX, bc_rank, &
      mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine 

  
  !> Broadcasts a character array from process root to all other processes in the group
  subroutine mpi_bcast_character( mpi_env, buffer, bcasting_rank )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer
    character(len=*), intent(inout) :: buffer(..)
    !> Broadcasting rank
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank
    bc_rank = mpi_env%root
    if ( present ( bcasting_rank ) ) bc_rank = bcasting_rank
    call mpi_bcast( buffer, len( buffer ) * size( buffer ), MPI_CHARACTER, bc_rank, &
      mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine 
end module
