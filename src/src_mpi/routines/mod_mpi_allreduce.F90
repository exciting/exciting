!> exciting wrappers for mpi_allreduce
module mod_mpi_allreduce
  use iso_c_binding, only: c_f_pointer, c_loc
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  ! A similar remark is found in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
  use mpi_f08, only: mpi_allreduce, mpi_in_place, mpi_sum, mpi_comm, MPI_INTEGER, &
    MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_Datatype, MPI_LOGICAL, MPI_LOR
#endif
  use precision, only: dp, i32, long_int

  implicit none
  private

  integer(i32), parameter :: big_int = huge(0_i32)/2

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
    logical, intent(inout) :: buffer(..)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_LOGICAL, MPI_LOR, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

! Generating mpi_allreduce_integer_i32( buffer, mpi_env, large_buffer_interface )
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#include "mpi_allreduce_template.inc"

! Generating mpi_allreduce_real_dp( buffer, mpi_env, large_buffer_interface )
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#include "mpi_allreduce_template.inc"

! Generating mpi_allreduce_complex_dp( buffer, mpi_env, large_buffer_interface )
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#include "mpi_allreduce_template.inc"

end module