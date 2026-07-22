!> exciting wrappers for mpi_bcast
module mod_mpi_bcast
  use mod_mpi_env, only: mpiinfo
  use iso_c_binding, only: c_loc, c_f_pointer ! Required for 1D mapping in < MPI 4.0
#ifdef MPI
  ! mpi_08 must be used to support Fortran 2008 and later
  use mpi_f08, only: mpi_bcast, mpi_comm, MPI_Datatype, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
    MPI_DOUBLE_COMPLEX, MPI_CHARACTER, MPI_LOGICAL
    
#ifdef USE_MPI4_LARGE_COUNTS
  use mpi_f08, only: MPI_COUNT_KIND
#endif

#endif
  use precision, only: dp, i32, long_int

  implicit none
  private

  integer(i32), parameter :: big_int = huge( 0_i32 ) / 2

  !> exciting wrapper of mpi_bcast 
  public :: xmpi_bcast

  interface xmpi_bcast
      module procedure :: &
        mpi_bcast_logical, &
        mpi_bcast_character, &
        mpi_bcast_integer_i32, &
        mpi_bcast_real_dp, &
        mpi_bcast_complex_dp
  end interface

contains

  !> Broadcasts logical array from the broadcasting_rank (root by default)
  subroutine mpi_bcast_logical( mpi_env, buffer, bcasting_rank )
    type(mpiinfo), intent(in) :: mpi_env
    logical, intent(inout) :: buffer(..)
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank, ierr
    
    bc_rank = mpi_env%root
    if ( present( bcasting_rank ) ) bc_rank = bcasting_rank

#ifdef USE_MPI4_LARGE_COUNTS
    call mpi_bcast( buffer, int( size( buffer, kind = long_int ), kind = MPI_COUNT_KIND), &
      MPI_LOGICAL, bc_rank, mpi_comm( mpi_env%comm ), ierr )
#else
    call mpi_bcast( buffer, size( buffer ), MPI_LOGICAL, bc_rank, &
      mpi_comm( mpi_env%comm ), ierr )
#endif

#endif
  end subroutine

  !> Broadcasts a character array from the broadcasting_rank (root by default)
  subroutine mpi_bcast_character( mpi_env, buffer, bcasting_rank )
    type(mpiinfo), intent(in) :: mpi_env
    character(len=*), intent(inout) :: buffer(..)
    integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
    integer(i32) :: bc_rank, ierr
    integer(long_int) :: total_chars
    
    bc_rank = mpi_env%root
    if ( present( bcasting_rank ) ) bc_rank = bcasting_rank
    
    total_chars = len( buffer, kind = long_int ) * size( buffer, kind = long_int )

#ifdef USE_MPI4_LARGE_COUNTS
    call mpi_bcast( buffer, int( total_chars, kind = MPI_COUNT_KIND ), MPI_CHARACTER, &
      bc_rank, mpi_comm( mpi_env%comm ), ierr )
#else
    call mpi_bcast( buffer, int( total_chars, kind = i32 ), MPI_CHARACTER, bc_rank, &
      mpi_comm( mpi_env%comm ), ierr )
#endif

#endif
  end subroutine 


! Generating mpi_bcast_integer_i32( mpi_env, buffer, bcasting_rank )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_bcast_template.inc"

! Generating mpi_bcast_real_dp( mpi_env, buffer, bcasting_rank )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_bcast_template.inc"

! Generating mpi_bcast_complex_dp( mpi_env, buffer, bcasting_rank )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_bcast_template.inc"

end module