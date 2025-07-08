!> exciting wrappers for `MPI_Reduce`
module mod_mpi_reduce
  use iso_c_binding, only: c_f_pointer, c_loc
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  ! A similar remark is found in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
  use mpi_f08, only: MPI_Reduce, mpi_in_place, mpi_sum, mpi_comm, MPI_INTEGER, &
    MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_Datatype
#endif
  use precision, only: dp, i32, long_int

  implicit none
  private

  integer(i32), parameter :: big_int = huge(0_i32)/2

  !> exciting wrapper of `MPI_Reduce`
  public :: xmpi_reduce

  interface xmpi_reduce
    module procedure :: mpi_reduce_integer_i32
    module procedure :: mpi_reduce_real_dp
    module procedure :: mpi_reduce_complex_dp
  end interface
contains

! Generating mpi_reduce_integer_i32( buffer, mpi_env, large_buffer_interface )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_reduce_template.inc"

! Generating mpi_reduce_real_dp( buffer, mpi_env, large_buffer_interface )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_reduce_template.inc"

! Generating mpi_reduce_complex_dp( buffer, mpi_env, large_buffer_interface )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_reduce_template.inc"
  
 
end module