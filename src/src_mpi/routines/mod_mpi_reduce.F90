!> exciting wrappers for `MPI_Reduce`
module mod_mpi_reduce
  use mod_mpi_env, only: mpiinfo
  use iso_c_binding, only: c_loc, c_f_pointer ! Required for 1D mapping in < MPI 4.0
#ifdef MPI
  ! mpi_08 must be used to support Fortran 2008 and later
  use mpi_f08, only: MPI_Reduce, MPI_IN_PLACE, mpi_sum, mpi_comm, MPI_INTEGER, &
    MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_Datatype
    
#ifdef USE_MPI4_LARGE_COUNTS
  use mpi_f08, only: MPI_COUNT_KIND
#endif

#endif
  use precision, only: dp, i32, long_int

  implicit none
  private

  integer(i32), parameter :: big_int = huge( 0_i32 ) / 2

  !> exciting wrapper of `MPI_Reduce`
  public :: xmpi_reduce

  interface xmpi_reduce
    module procedure :: mpi_reduce_integer_i32
    module procedure :: mpi_reduce_real_dp
    module procedure :: mpi_reduce_complex_dp
  end interface
contains

! Generating mpi_reduce_integer_i32( buffer, mpi_env )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_reduce_template.inc"

! Generating mpi_reduce_real_dp( buffer, mpi_env )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_reduce_template.inc"

! Generating mpi_reduce_complex_dp( buffer, mpi_env )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_reduce_template.inc"
  
end module