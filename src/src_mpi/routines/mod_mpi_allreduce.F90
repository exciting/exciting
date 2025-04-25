!> exciting wrappers for mpi_allreduce
module mod_mpi_allreduce
#ifdef MPI
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  ! A similar remark is found in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
  use mpi_f08, only: mpi_allreduce, mpi_in_place, mpi_sum, mpi_comm, MPI_INTEGER, &
    MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX
#endif
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

  subroutine mpi_allreduce_integer_dp( buffer, mpi_env )
    integer(i32), intent(inout) :: buffer(..)
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_INTEGER, mpi_sum, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

  subroutine mpi_allreduce_real_dp( buffer, mpi_env )
    real(dp), intent(inout) :: buffer(..)
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

  subroutine mpi_allreduce_complex_dp( buffer, mpi_env )
    complex(dp), intent(inout) :: buffer(..)
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_DOUBLE_COMPLEX, mpi_sum, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

end module