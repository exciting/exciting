!> exciting wrappers for mpi_allreduce
module mod_mpi_allreduce
#ifdef MPI
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  ! A similar remark is found in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
  use mpi_f08, only: mpi_allreduce, mpi_in_place, mpi_sum, mpi_comm, MPI_INTEGER, &
    MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_LOGICAL, MPI_LOR
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

  !> Wrapper for mpi_allreduce targeting `integer(i32)` data. 
  !> Sum values from all processes and distribute the result back to all processes.
  subroutine mpi_allreduce_integer_dp( buffer, mpi_env )
    !> Data
    integer(i32), intent(inout) :: buffer(..)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_INTEGER, mpi_sum, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

  !> Wrapper for mpi_allreduce targeting `real(dp)` data.
  !> Sum values from all processes and distribute the result back to all processes.
  subroutine mpi_allreduce_real_dp( buffer, mpi_env )
    !> Data
    real(dp), intent(inout) :: buffer(..)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_DOUBLE_PRECISION, mpi_sum, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

  !> Wrapper for mpi_allreduce targeting `complex(dp)` data.
  !> Sum values from all processes and distribute the result back to all processes.
  subroutine mpi_allreduce_complex_dp( buffer, mpi_env )
    !> Data
    complex(dp), intent(inout) :: buffer(..)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
    integer(i32) :: ierr
    
    call mpi_allreduce( mpi_in_place, buffer, size( buffer ), MPI_DOUBLE_COMPLEX, mpi_sum, mpi_comm(mpi_env%comm), ierr )
#endif
  end subroutine

end module