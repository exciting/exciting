!> Module for mpi_allgather wrappers.
module mod_mpi_allgather
  use mod_mpi_env, only: mpiinfo
  use precision, only: dp, i32

#ifdef MPI
  use mpi_f08, only: MPI_INTEGER, MPI_DOUBLE_COMPLEX, mpi_allgather, mpi_allgatherv, &
    mpi_in_place, MPI_DATATYPE_NULL, mpi_comm, MPI_DOUBLE_PRECISION
#endif

  implicit none 
  private
  
  public :: xmpi_allgather, xmpi_allgatherv

  !> Wrappers for mpi_allgether.
  interface xmpi_allgather
    module procedure :: &
      mpi_allgather_integer_i32
  end interface

  !> Wrappers for mpi_allgetherv.
  interface xmpi_allgatherv
    module procedure :: &
      mpi_allgatherv_in_place_complex_dp, &
      mpi_allgatherv_in_place_real_dp
  end interface

contains

  !> Wrapper for mpi_allgather for an `integer(i32)` scalar.
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgather_integer_i32( mpi_env, send_buffer, receive_buffer )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Integer to be send by the current rank.
    integer(i32), intent(in) :: send_buffer
    !> Array that contains the the integers, recieved by all ranks
    integer(i32), allocatable, intent(out) :: receive_buffer(:)
#ifdef MPI    

    allocate( receive_buffer( mpi_env%procs ) )
    receive_buffer( mpi_env%rank + 1) = send_buffer
    call mpi_allgather( send_buffer, 1, MPI_INTEGER, receive_buffer, 1, MPI_INTEGER, &
      mpi_comm( mpi_env%comm ), mpi_env%ierr )
#else
    allocate( receive_buffer(1), source=send_buffer )
#endif
  end subroutine mpi_allgather_integer_i32


  !> Wrapper for mpi_allgatherv for a `complex(dp)` array.
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgatherv_in_place_complex_dp( mpi_env, buffer, chunk_size )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer. On input, contains the data of the current rank
    !> On output, contains the data of all ranks.
    complex(dp), intent(inout) :: buffer(..)
    !> Number of elements handled by the current rank.
    integer(i32), intent(in) :: chunk_size
#ifdef MPI
    integer(i32), allocatable :: receive_counts(:), displacements(:)

    call xmpi_allgather( mpi_env, chunk_size, receive_counts )
    call calculate_displacements( mpi_env, receive_counts, displacements )
    
    call mpi_allgatherv( mpi_in_place, 0, MPI_DATATYPE_NULL, buffer, receive_counts, &
      displacements, MPI_DOUBLE_COMPLEX, mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine mpi_allgatherv_in_place_complex_dp

  !> Wrapper for mpi_allgatherv for a `real(dp)` array.
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgatherv_in_place_real_dp( mpi_env, buffer, chunk_size )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer. On input, contains the data of the current rank
    !> On output, contains the data of all ranks.
    real(dp), intent(inout) :: buffer(..)
    !> Number of elements handled by the current rank.
    integer(i32), intent(in) :: chunk_size
#ifdef MPI
    integer(i32), allocatable :: receive_counts(:), displacements(:)

    call xmpi_allgather( mpi_env, chunk_size, receive_counts )
    call calculate_displacements( mpi_env, receive_counts, displacements )
    
    call mpi_allgatherv( mpi_in_place, 0, MPI_DATATYPE_NULL, buffer, receive_counts, &
      displacements, MPI_DOUBLE_PRECISION, mpi_comm( mpi_env%comm ), mpi_env%ierr )
#endif
  end subroutine mpi_allgatherv_in_place_real_dp


  !> Calculate the number of elements, the data chunk of each rank is displaced by.
  subroutine calculate_displacements( mpi_env, receive_counts, displacements )
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Array that holds the size of the data chunk for each rank.
    integer(i32), intent(in) :: receive_counts(:)
    !> Number of elements, the data chunk of each rank is displaced by.
    integer(i32), allocatable, intent(out) :: displacements(:)

    integer(i32) :: rank

    allocate( displacements(mpi_env%procs) )
    do rank = 0, mpi_env%procs - 1
      displacements(rank + 1) = sum( receive_counts(1 : rank) )
    end do
  end subroutine 
  

end module mod_mpi_allgather