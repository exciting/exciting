!> Serial routine overloads for exciting's MPI wrappers.
!> These typically do nothing other than ensure the code
!> will compile in serial, without needing to dress all MPI 
!> calls throughout the code in preprocessor variables. 
module mod_serial_bcast
    use mod_mpi_env, only: mpiinfo
    use precision, only: sp, dp

    implicit none
    private

    !> exciting wrapper of mpi_bcast
    public :: xmpi_bcast

    interface xmpi_bcast
      module procedure :: mpi_bcast_rank0_int_sp, &
                & mpi_bcast_rank0_logical, &
                & mpi_bcast_rank0_real_dp, mpi_bcast_rank1_real_dp, mpi_bcast_rank2_real_dp, mpi_bcast_rank3_real_dp, &
                & mpi_bcast_rank0_complex_dp, mpi_bcast_rank2_complex_dp, mpi_bcast_rank3_complex_dp, mpi_bcast_rank4_complex_dp, &
                & mpi_bcast_character, mpi_bcast_character_array 
    end interface 

contains    

    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank0_logical(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      logical,  intent(in) :: buffer
    end subroutine 

    
    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank0_int_sp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env 
      !> Buffer
      integer(sp), intent(in) :: buffer
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank0_real_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env 
      !> Buffer
      real(dp), intent(in) :: buffer
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank1_real_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env 
      !> Buffer
      real(dp), intent(in) :: buffer(:)
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank2_real_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      real(dp), intent(in) :: buffer(:,:) 
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank3_real_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      real(dp), intent(in) :: buffer(:,:,:) 
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank0_complex_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      complex(dp), intent(in) :: buffer
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank1_complex_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      complex(dp), intent(in) :: buffer(:)
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank2_complex_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      complex(dp), intent(in) :: buffer(:,:)
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank3_complex_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      complex(dp), intent(in) :: buffer(:,:,:)
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_rank4_complex_dp(mpi_env, buffer)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env
      !> Buffer
      complex(dp), intent(in) :: buffer(:,:,:,:)
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_character(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env 
        !> Buffer 
        character(len=*), intent(in):: buffer
    end subroutine 


    !> Dummy routine for serial version of mpi_bcast
    subroutine mpi_bcast_character_array(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env 
        !> Buffer 
        character(len=*), intent(in):: buffer(:)
    end subroutine 

end module  
