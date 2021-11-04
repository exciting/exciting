!> exciting wrappers for mpi_bcast
module mod_mpi_bcast   
#ifdef MPI
      use mpi
      use mod_mpi_env, only: mpiinfo
      use precision, only: sp, dp 
  
      implicit none
      private
  
      !> exciting wrapper of mpi_bcast 
      public :: xmpi_bcast
  
      interface xmpi_bcast
          module procedure :: mpi_bcast_rank0_int_sp, mpi_bcast_rank0_real_dp,&
            mpi_bcast_rank1_real_dp, mpi_bcast_rank2_real_dp,&
            & mpi_bcast_rank0_complex_dp, mpi_bcast_rank2_complex_dp, &
            & mpi_bcast_rank3_complex_dp, mpi_bcast_rank4_complex_dp, &
            & mpi_bcast_character, mpi_bcast_character_array 
      end interface 
  
  contains
  
      !> Broadcasts integer scalar from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank0_int_sp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        integer(sp),  intent(in) :: buffer
        call MPI_BCAST(buffer, 1, MPI_INTEGER, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts real scalar from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank0_real_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        real(dp), intent(in) :: buffer
  
        call MPI_BCAST(buffer,1, MPI_DOUBLE_PRECISION, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts real array (Rank = 1) from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank1_real_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        real(dp), intent(in) :: buffer(:) 
  
        call MPI_BCAST(buffer, size(buffer), MPI_DOUBLE_PRECISION, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts real array (Rank = 2) from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank2_real_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        real(dp), intent(in) :: buffer(:,:) 
  
        call MPI_BCAST(buffer, size(buffer), MPI_DOUBLE_PRECISION, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts complex scalar from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank0_complex_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        complex(dp), intent(in) :: buffer
  
        call MPI_BCAST(buffer, 1, MPI_DOUBLE_COMPLEX, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
      !> Broadcasts complex array (Rank = 1) from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank1_complex_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        complex(dp), intent(in) :: buffer(:)
  
        call MPI_BCAST(buffer, size(buffer), MPI_DOUBLE_COMPLEX, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts complex array (Rank = 2) from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank2_complex_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        complex(dp), intent(in) :: buffer(:,:)
  
        call MPI_BCAST(buffer, size(buffer), MPI_DOUBLE_COMPLEX, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts complex array (Rank = 3) from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank3_complex_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        complex(dp), intent(in) :: buffer(:,:,:)
  
        call MPI_BCAST(buffer, size(buffer), MPI_DOUBLE_COMPLEX, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
      !> Broadcasts complex array (Rank = 4) from the process with rank root to all other processes of the group. 
      subroutine mpi_bcast_rank4_complex_dp(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Buffer
        complex(dp), intent(in) :: buffer(:,:,:,:)
  
        call MPI_BCAST(buffer, size(buffer), MPI_DOUBLE_COMPLEX, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts character from process root to all other processes in the group
      subroutine mpi_bcast_character(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Character buffer
        character(len=*), intent(in) :: buffer 
        call MPI_BCAST(buffer, len(buffer), MPI_CHARACTER, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
  
  
      !> Broadcasts a character array (Rank = 1) from process root to all other processes in the group
      subroutine mpi_bcast_character_array(mpi_env, buffer)
        !> MPI environment
        type(mpiinfo), intent(inout) :: mpi_env
        !> Character buffer
        character(len=*), intent(in) :: buffer(:) 
        call MPI_BCAST(buffer, len(buffer) * size(buffer), MPI_CHARACTER, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      end subroutine 
#endif
end module
