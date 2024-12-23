module mod_serial_gather
#ifndef MPI
    use mod_mpi_env, only: mpiinfo
    use precision, only: dp, i32
  
    implicit none 
  
    private
    
    public :: xmpi_gather, xmpi_gatherv
  
    !> exciting wrapper for `mpi_gather`
    interface xmpi_gather
      module procedure :: &
        mpi_gather_rank1_int_sp
    end interface
  
    !> exciting wrapper for `mpi_gatherv`
    interface xmpi_gatherv
      module procedure :: &
        mpi_gatherv_rank1_int_sp, &
        mpi_gatherv_rank2_real_dp, &
        mpi_gatherv_rank3_complex_dp
    end interface
  
  contains
  
    !> Gather int arrays of rank 1. The arrays must have same sizes
    subroutine mpi_gather_rank1_int_sp( mpi_env, array_to_send, receive_buffer )
      !> MPI environment
      type(mpiinfo), intent(in) :: mpi_env
      !> Array to send to the root rank
      integer(i32), contiguous, intent(in) :: array_to_send(:)
      !> Array where all sent arrays are gathered (relevant only for the root rank)
      integer(i32), allocatable, intent(out) :: receive_buffer(:)

      receive_buffer = array_to_send
    end subroutine
  
    !> Gather int arrays of rank 1. The arrays may have different sizes
    subroutine mpi_gatherv_rank1_int_sp( mpi_env, array_to_send, receive_buffer )
      !> MPI environment
      type(mpiinfo), intent(in) :: mpi_env
      !> Array to send to the root rank
      integer(i32), intent(in) :: array_to_send(:)
      !> Array where all sent arrays are gathered (relevant only for the root rank)
      integer(i32), allocatable, intent(out) :: receive_buffer(:)

      receive_buffer = array_to_send
    end subroutine
  
    
    !> Gather real-dp array of rank 2. The arrays may have different sizes
    subroutine mpi_gatherv_rank2_real_dp( mpi_env, array_to_send, receive_buffer )
      !> MPI environment
      type(mpiinfo), intent(in) :: mpi_env
      !> Array to send to the root rank
      real(dp), intent(in) :: array_to_send(:, :)
      !> Array where all sent arrays are gathered (relevant only for the root rank)
      real(dp), allocatable, intent(out) :: receive_buffer(:, :)

      receive_buffer = array_to_send
    end subroutine
    

    !> Gather complex-dp array of rank 3. The arrays may have different sizes
    subroutine mpi_gatherv_rank3_complex_dp( mpi_env, array_to_send, receive_buffer )
      !> MPI environment
      type(mpiinfo), intent(in) :: mpi_env
      !> Array to send to the root rank
      complex(dp), intent(in) :: array_to_send(:, :, :)
      !> Array where all sent arrays are gathered (relevant only for the root rank)
      complex(dp), allocatable, intent(out) :: receive_buffer(:, :, :)

      receive_buffer = array_to_send
    end subroutine
#endif   
end module