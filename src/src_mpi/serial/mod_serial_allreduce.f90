!> Serial routine overloads for exciting's MPI wrappers.
!> These typically do nothing other than ensure the code
!> will compile in serial, without needing to dress all MPI 
!> calls throughout the code in preprocessor variables.
module mod_serial_allreduce
#ifndef MPI
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
end subroutine

subroutine mpi_allreduce_real_dp( data, mpi_env )
  real(dp), intent(inout) :: data(..)
  type(mpiinfo), intent(in) :: mpi_env
end subroutine

subroutine mpi_allreduce_complex_dp( data, mpi_env )
  complex(dp), intent(inout) :: data(..)
  type(mpiinfo), intent(in) :: mpi_env
end subroutine
#endif
end module