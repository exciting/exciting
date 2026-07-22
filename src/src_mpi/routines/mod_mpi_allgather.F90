!> Module for mpi_allgather(v) wrappers.
module mod_mpi_allgather
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  ! mpi_08 must be used to support Fortran 2008 and later
  use mpi_f08, only: mpi_allgather, mpi_allgatherv, MPI_IN_PLACE, mpi_comm, &
    MPI_Datatype, MPI_DATATYPE_NULL, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, &
    MPI_INTEGER, MPI_REAL, MPI_INTEGER8
  use modmpi, only: terminate_mpi_env
#ifdef USE_MPI4_LARGE_COUNTS
  use mpi_f08, only: MPI_COUNT_KIND, MPI_ADDRESS_KIND
#endif

#endif
  use precision, only: dp, sp, i32, long_int

  implicit none
  private
  
  public :: xmpi_allgather, xmpi_allgatherv

  !> Wrappers for mpi_allgather.
  interface xmpi_allgather
    module procedure :: &
      mpi_allgather_integer_i32, &
      mpi_allgather_integer_long_int
  end interface

  !> Wrappers for mpi_allgatherv.
  interface xmpi_allgatherv
    module procedure :: &
      mpi_allgatherv_in_place_complex_dp, &
      mpi_allgatherv_in_place_real_dp, &
      mpi_allgatherv_in_place_real_sp, &
      mpi_allgatherv_in_place_integer_i32, &
      mpi_allgatherv_in_place_complex_dp_large_buffer, &
      mpi_allgatherv_in_place_real_dp_large_buffer, &
      mpi_allgatherv_in_place_real_sp_large_buffer, &
      mpi_allgatherv_in_place_integer_i32_large_buffer
  end interface

contains

  !> Wrapper for mpi_allgather for an `integer(long_int)` scalar.
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgather_integer_long_int( mpi_env, send_buffer, receive_buffer )
    type(mpiinfo), intent(in) :: mpi_env
    integer(long_int), intent(in) :: send_buffer
    integer(long_int), allocatable, intent(out) :: receive_buffer(:)
#ifdef MPI
    integer(i32) :: ierr
    allocate( receive_buffer( mpi_env%procs ) )
    call mpi_allgather( send_buffer, 1, MPI_INTEGER8, receive_buffer, &
      1, MPI_INTEGER8, mpi_comm( mpi_env%comm ), ierr )
#else
    allocate( receive_buffer(1), source=send_buffer )
#endif
  end subroutine mpi_allgather_integer_long_int

  !> Wrapper for mpi_allgather for an `integer(i32)` scalar.
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgather_integer_i32( mpi_env, send_buffer, receive_buffer )
    type(mpiinfo), intent(in) :: mpi_env
    integer(i32), intent(in) :: send_buffer
    integer(i32), allocatable, intent(out) :: receive_buffer(:)
#ifdef MPI
    integer(i32) :: ierr
    allocate( receive_buffer( mpi_env%procs ) )
    call mpi_allgather( send_buffer, 1, MPI_INTEGER, receive_buffer, 1, MPI_INTEGER, &
      mpi_comm( mpi_env%comm ), ierr )
#else
    allocate( receive_buffer(1), source=send_buffer )
#endif
  end subroutine mpi_allgather_integer_i32

! Generating mpi_allgatherv_integer_i32( mpi_env, buffer, chunk_size )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_allgatherv_template.inc"

! Generating mpi_allgatherv_real_dp( mpi_env, buffer, chunk_size )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_allgatherv_template.inc"

! Generating mpi_allgatherv_real_sp( mpi_env, buffer, chunk_size )
#define TYPE1 real
#define PRECISION1 sp
#define MPIDATATYPE1 MPI_REAL
#include "mpi_allgatherv_template.inc"

! Generating mpi_allgatherv_complex_dp( mpi_env, buffer, chunk_size )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_allgatherv_template.inc"

end module mod_mpi_allgather