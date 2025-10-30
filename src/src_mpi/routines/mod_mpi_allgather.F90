!> Module for mpi_allgatherv wrappers.
module mod_mpi_allgather
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  use mpi_f08, only: mpi_allgatherv, mpi_in_place, mpi_comm, MPI_DATATYPE_NULL, &
    MPI_Datatype, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_INTEGER, MPI_REAL 
#endif
  use precision, only: i32, dp, sp

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
      mpi_allgatherv_in_place_real_dp, &
      mpi_allgatherv_in_place_real_sp, &
      mpi_allgatherv_in_place_integer_i32
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

! Generating mpi_allgatherv_integer_i32( mpi_env, buffer, chunk_size )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_allgatherv_template.inc"

! Generating mpi_allreduce_real_dp( mpi_env, buffer, chunk_size )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_allgatherv_template.inc"

! Generating mpi_allreduce_real_sp( mpi_env, buffer, chunk_size )
#define TYPE1 real
#define PRECISION1 sp
#define MPIDATATYPE1 MPI_REAL
#include "mpi_allgatherv_template.inc"

! Generating mpi_allreduce_complex_dp( mpi_env, buffer, chunk_size )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_allgatherv_template.inc"

end module mod_mpi_allgather