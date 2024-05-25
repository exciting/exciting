!> exciting wrappers for mpi_allreduce
!> Every addition to a parallel wrapper also requires a serial overload
module mod_mpi_allreduce
#ifdef MPI
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  ! A similar remark is found in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
  use mpi_f08, only: mpi_allreduce, MPI_IN_PLACE, MPI_SUM, MPI_COMM, &
                     MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX

  use mod_mpi_env, only: mpiinfo
  use precision, only: dp, i32

  implicit none

  private

  !> exciting wrapper of mpi_allreduce
  public :: xmpi_allreduce

  interface xmpi_allreduce
    module procedure :: mpi_allreduce_integer_dp
    module procedure :: mpi_allreduce_real_dp
    module procedure :: mpi_allreduce_complex_dp
  end interface
contains

subroutine mpi_allreduce_integer_dp( data, mpi_env )
  integer(i32), intent(inout) :: data(..)
  type(mpiinfo), intent(in) :: mpi_env

  integer(i32) :: ierr
  
  call mpi_allreduce( MPI_IN_PLACE, data, size( data ), MPI_INTEGER, MPI_SUM, MPI_COMM(mpi_env%comm), ierr )
end subroutine

subroutine mpi_allreduce_real_dp( data, mpi_env )
  real(dp), intent(inout) :: data(..)
  type(mpiinfo), intent(in) :: mpi_env

  integer(i32) :: ierr
  
  call mpi_allreduce( MPI_IN_PLACE, data, size( data ), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM(mpi_env%comm), ierr )
end subroutine

subroutine mpi_allreduce_complex_dp( data, mpi_env )
  complex(dp), intent(inout) :: data(..)
  type(mpiinfo), intent(in) :: mpi_env

  integer(i32) :: ierr
  
  call mpi_allreduce( MPI_IN_PLACE, data, size( data ), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM(mpi_env%comm), ierr )
end subroutine

#endif
end module