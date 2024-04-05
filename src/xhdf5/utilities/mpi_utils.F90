!> MPI utilities for the HDF5 wrappers
module mpi_utils
  
#ifdef MPI  
  use mpi
#endif

  use trace, only: trace_back
  use iso_fortran_env, only: error_unit

  
  implicit none

 
  !> Rank of the root process.
  integer, parameter, public :: root_rank = 0

  !> Type for handling the MPI communicator.
  type mpi_comm_type
    integer :: handle
  end type 



  contains 

  !> Return the rank an MPI communicator handle belongs to.
  !> 
  !> Wrapper for `mpi_comm_rank`.
  integer function comm_to_rank(mpi_comm)
    !> MPI communicator handle
    type(mpi_comm_type), intent(in) :: mpi_comm

    integer :: mpierr

#ifdef MPI
    call mpi_comm_rank(mpi_comm%handle, comm_to_rank, mpierr)
    call handle_mpi_error(mpi_comm, 'mpi_comm_rank', mpierr)
#else 
    comm_to_rank = root_rank
#endif
    
  end function comm_to_rank


  !> Terminate an MPI environment.
  subroutine terminate_all_mpi_processes(mpi_comm, message)
    !> MPI communicator.
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Error message
    character(len=*), optional, intent(in) :: message

    !> Error code to send to `mpi_abort`.
    integer, parameter :: error_code = 101

    integer :: mpierr

#ifdef MPI
    if(comm_to_rank(mpi_comm) == root_rank) then
      if(present(message)) write(error_unit, *) trim(adjustl(message))
      call trace_back()
    end if
    call mpi_barrier(MPI_COMM_WORLD, mpierr)
    call mpi_abort(MPI_COMM_WORLD, error_code, mpierr)
#else
    if(present(message)) write(error_unit, *) trim(adjustl(message))
    call trace_back()
    stop
#endif

  end subroutine terminate_all_mpi_processes


  !> Assert that the MPI error flag `mpierr` is zero. 
  !> When not, `terminate_mpi_comm` is called to terminate the MPI environment.
  subroutine handle_mpi_error(mpi_comm, calling_routine, mpierr)
    !> MPI communicator.
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Name of the hdf5 routine.
    character(*), intent(in) :: calling_routine
    !> MPI error flag.
    integer, intent(in) :: mpierr

    character(3) :: mpierr_char
    character(:), allocatable :: message

    if (mpierr /= 0) then
      write(mpierr_char, '(I3)') mpierr
      message = 'MPI routine ' // calling_routine // ' failed with mpierr = '// mpierr_char
      call terminate_all_mpi_processes(mpi_comm, message)
    end if

  end subroutine handle_mpi_error

end module 