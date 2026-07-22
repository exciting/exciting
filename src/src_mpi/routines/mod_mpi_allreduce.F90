!> exciting wrappers for mpi_allreduce
module mod_mpi_allreduce
  use mod_mpi_env, only: mpiinfo
  use iso_c_binding, only: c_loc, c_f_pointer ! Required for 1D mapping in < MPI 4.0
#ifdef MPI
  ! mpi_08 must be used to support Fortran 2008 and later
  use mpi_f08, only: mpi_allreduce, MPI_IN_PLACE, mpi_sum, mpi_comm, MPI_INTEGER, &
    MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_Datatype, MPI_LOGICAL, MPI_LOR
#ifdef USE_MPI4_LARGE_COUNTS
  use mpi_f08, only: MPI_COUNT_KIND
#endif

#endif
  use precision, only: dp, i32, long_int

  implicit none
  private

  integer(long_int), parameter :: big_int = int(huge(0_i32), kind = long_int) / 2_long_int

  !> exciting wrapper of mpi_allreduce
  public :: xmpi_allreduce

  interface xmpi_allreduce
    module procedure :: mpi_allreduce_integer_i32
    module procedure :: mpi_allreduce_real_dp
    module procedure :: mpi_allreduce_complex_dp
    module procedure :: mpi_allreduce_logical
  end interface
contains
  
  !> Wrapper for mpi_allreduce targeting `logical` data.
  !> Sum (logical OR) values from all processes and distribute the result back to all processes.
  subroutine mpi_allreduce_logical( buffer, mpi_env )
    !> Data
    logical, contiguous, target, intent(inout) :: buffer(..)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(long_int) :: buffer_size
    integer(long_int) :: first, last, chunk_elements
    logical, pointer :: buf_1d(:)
    integer(i32) :: ierr
    
    buffer_size = size( buffer, kind = long_int )
    
    if ( buffer_size > big_int ) then
#ifdef USE_MPI4_LARGE_COUNTS
      call mpi_allreduce( MPI_IN_PLACE, buffer, int(buffer_size, kind = MPI_COUNT_KIND), &
                          MPI_LOGICAL, MPI_LOR, mpi_comm(mpi_env%comm), ierr )
#else
      call c_f_pointer( c_loc(buffer), buf_1d, [buffer_size] )
      
      first = 1_long_int
      do while ( first <= buffer_size )
        last = min( first + big_int - 1_long_int, buffer_size )
        chunk_elements = last - first + 1_long_int
        
        call mpi_allreduce( MPI_IN_PLACE, buf_1d(first:last), int( chunk_elements, kind = i32 ), &
                            MPI_LOGICAL, MPI_LOR, mpi_comm(mpi_env%comm), ierr )
        first = last + 1_long_int
      end do
#endif
    else
      call mpi_allreduce( MPI_IN_PLACE, buffer, int( buffer_size, kind = i32 ), &
                          MPI_LOGICAL, MPI_LOR, mpi_comm(mpi_env%comm), ierr )
    end if
#endif
  end subroutine

! Generating mpi_allreduce_integer_i32( buffer, mpi_env )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_allreduce_template.inc"

! Generating mpi_allreduce_real_dp( buffer, mpi_env )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_allreduce_template.inc"

! Generating mpi_allreduce_complex_dp( buffer, mpi_env )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_allreduce_template.inc"

end module