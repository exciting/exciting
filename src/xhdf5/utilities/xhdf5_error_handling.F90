!> Error handling for xhdf5.
module xhdf5_error_handling

  use mpi_utils, only: terminate_all_mpi_processes, mpi_comm_type

  implicit none

  private
  public :: xhdf5_assert, handle_hdf5_error, abort_if_not_hdf5


  contains
  

  !> Assert if a condition is true. If not the all processes of the MPI comminicator will be terminated.
  subroutine xhdf5_assert(mpi_comm, condition, message)
    !> MPI communicator.
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Condition that must be true to pass the xhdf5_assert.
    logical, intent(in) :: condition
    !> Message to print if assertion fails
    character(*), intent(in) :: message

    if (.not. condition) then
        call terminate_all_mpi_processes(mpi_comm, message=message)
    end if

  end subroutine xhdf5_assert


  !> Assert that the HDF5 error flag `h5err` is zero. 
  !> If not, call [[terminate_mpi_comm]] with an error message containing the name of `calling_routine`.
  subroutine handle_hdf5_error(mpi_comm, calling_routine, h5err)
    !> MPI communicator.
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Name of the hdf5 routine.
    character(*), intent(in) :: calling_routine
    !> HDF5 error flag.
    integer, intent(in) :: h5err

    character(3) :: h5err_char
    character(:), allocatable :: message

    if (h5err /= 0) then
      write(h5err_char, '(I3)') h5err
      message = 'HDF5 routine ' // calling_routine // ' failed with h5err = '// h5err_char
      call terminate_all_mpi_processes(mpi_comm, message)
    end if
  end subroutine handle_hdf5_error
  

  !> Abort run if exciting was compiled exciting without HDF5.
  !> Call this routine in features that require HDF5 in the beginning.
  ! This should go to an module in exciting
  subroutine abort_if_not_hdf5(mpi_env, message)
    use modmpi, only: mpiinfo, terminate_mpi_env, mpiglobal
    !> MPI environment to terminate.
    type(mpiinfo), intent(inout), optional :: mpi_env
    !> Message to print to the terminal
    character(*), intent(in), optional :: message

    character(:), allocatable :: message_local
    type(mpiinfo) :: mpi_env_local

    character(*), parameter :: default_message = "Error: exciting needs to be linked to HDF5 to run this workflow."
    
    message_local = default_message
    if(present(message)) message_local = message

    mpi_env_local = mpiglobal
    if(present(mpi_env)) mpi_env_local = mpi_env

#ifndef _HDF5_ 
    call terminate_mpi_env(mpi_env_local, message=message)
#endif  
  end subroutine abort_if_not_hdf5

end module xhdf5_error_handling  