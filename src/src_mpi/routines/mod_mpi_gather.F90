module mod_mpi_gather
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
  ! assumed rank arrays in calls to openmpi-MPI subroutines.
  ! A similar remark is found in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
#ifdef MPI
  use mpi_f08, only: mpi_gather, mpi_gatherv, mpi_comm, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, MPI_INTEGER
#endif
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
#ifdef MPI
    integer(i32) :: fake_buffer(1), i_error

    if( mpi_env%is_root ) then
      allocate( receive_buffer((mpi_env%procs)*size( array_to_send )) )
      call mpi_gather( array_to_send, size( array_to_send ), MPI_INTEGER, receive_buffer, &
        size( array_to_send ), MPI_INTEGER, mpi_env%root, mpi_comm( mpi_env%comm ), i_error )
    else
      call mpi_gather( array_to_send, size( array_to_send ), MPI_INTEGER, fake_buffer, &
        size( array_to_send ), MPI_INTEGER, mpi_env%root, mpi_comm( mpi_env%comm ), i_error )
    end if
#else
    receive_buffer = array_to_send
#endif
  end subroutine

  !> Gather int arrays of rank 1. The arrays may have different sizes
  subroutine mpi_gatherv_rank1_int_sp( mpi_env, array_to_send, receive_buffer )
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Array to send to the root rank
    integer(i32), intent(in) :: array_to_send(:)
    !> Array where all sent arrays are gathered (relevant only for the root rank)
    integer(i32), allocatable, intent(out) :: receive_buffer(:)
#ifdef MPI    
    real(dp) :: fake_buffer(1)
    integer(i32) :: fake_counts(1), fake_displacements(1), i_error
    integer(i32), allocatable :: receive_counts(:), receive_displacements(:)

    call get_receive_counts( mpi_env, size( array_to_send ), receive_counts )
    if( .not. mpi_env%is_root ) then
      call mpi_gatherv( array_to_send, size( array_to_send ), MPI_INTEGER, &
        fake_buffer, fake_counts, fake_displacements, MPI_INTEGER, mpi_env%root, &
        mpi_comm( mpi_env%comm ), i_error )
    else
      call get_displacements( receive_counts, receive_displacements )
      allocate( receive_buffer( sum( receive_counts ) ) )
      call mpi_gatherv( array_to_send, size( array_to_send ), MPI_INTEGER, &
        receive_buffer, receive_counts, receive_displacements, MPI_INTEGER, mpi_env%root, &
        mpi_comm( mpi_env%comm ), i_error )
    end if
#else
    receive_buffer = array_to_send
#endif
  end subroutine

  
  !> Gather real-dp array of rank 2. The arrays may have different sizes
  subroutine mpi_gatherv_rank2_real_dp( mpi_env, array_to_send, receive_buffer )
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Array to send to the root rank
    real(dp), intent(in) :: array_to_send(:, :)
    !> Array where all sent arrays are gathered (relevant only for the root rank)
    real(dp), allocatable, intent(out) :: receive_buffer(:, :)
#ifdef MPI
    real(dp) :: fake_buffer(1)
    integer(i32) :: m, n, fake_counts(1), fake_displacements(1), i_error
    integer(i32), allocatable :: receive_counts(:), receive_displacements(:)

    call get_receive_counts( mpi_env, size( array_to_send ), receive_counts )
    if( .not. mpi_env%is_root ) then
      call mpi_gatherv( array_to_send, size( array_to_send ), MPI_DOUBLE_PRECISION, &
        fake_buffer, fake_counts, fake_displacements, MPI_DOUBLE_PRECISION, mpi_env%root, &
        mpi_comm( mpi_env%comm ), i_error )
    else
      call get_displacements( receive_counts, receive_displacements )
      m = size( array_to_send, 1 )
      n = int( sum( receive_counts ) / m, kind = i32 )
      allocate( receive_buffer(m, n) )
      call mpi_gatherv( array_to_send, size( array_to_send ), MPI_DOUBLE_PRECISION, &
        receive_buffer, receive_counts, receive_displacements, MPI_DOUBLE_PRECISION, &
        mpi_env%root, mpi_comm( mpi_env%comm ), i_error )
    end if
#else
    receive_buffer = array_to_send
#endif
  end subroutine

  
  !> Gather complex-dp array of rank 3. The arrays may have different sizes
  subroutine mpi_gatherv_rank3_complex_dp( mpi_env, array_to_send, receive_buffer )
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Array to send to the root rank
    complex(dp), intent(in) :: array_to_send(:, :, :)
    !> Array where all sent arrays are gathered (relevant only for the root rank)
    complex(dp), allocatable, intent(out) :: receive_buffer(:, :, :)
#ifdef MPI
    complex(dp) :: fake_buffer(1)
    integer(i32) :: m, n, k, fake_counts(1), fake_displacements(1), i_error
    integer(i32), allocatable :: receive_counts(:), receive_displacements(:)

    call get_receive_counts( mpi_env, size( array_to_send ), receive_counts )
    if( .not. mpi_env%is_root ) then
      call mpi_gatherv( array_to_send, size( array_to_send ), MPI_DOUBLE_COMPLEX, &
        fake_buffer, fake_counts, fake_displacements, MPI_DOUBLE_COMPLEX, mpi_env%root, &
        mpi_comm( mpi_env%comm ), i_error )
    else
      call get_displacements( receive_counts, receive_displacements )
      m = size( array_to_send, 1 )
      n = size( array_to_send, 2 )
      k = int( sum( receive_counts ) / (m * n), kind = i32 )
      allocate( receive_buffer(m, n, k) )
      call mpi_gatherv( array_to_send, size( array_to_send ), MPI_DOUBLE_COMPLEX, &
        receive_buffer, receive_counts, receive_displacements, MPI_DOUBLE_COMPLEX, &
        mpi_env%root, mpi_comm( mpi_env%comm ), i_error )
    end if
#else
    receive_buffer = array_to_send
#endif
  end subroutine


  !> (Private) Wrapper to [[xmpi_gather]] with a better name of what is expected to do
  subroutine get_receive_counts( mpi_env, send_count, receive_counts )
    type(mpiinfo), intent(in) :: mpi_env
    integer(i32), intent(in) :: send_count
    integer(i32), allocatable, intent(out) :: receive_counts(:)

    call xmpi_gather( mpi_env, [send_count], receive_counts )
  end subroutine

  !> (Private) Get displacements from counts
  subroutine get_displacements( receive_counts, displacements )
    !> Number of elements that each rank sents to the root rank
    integer(i32), intent(in) :: receive_counts(:)
    !> Displacement of each array, when gathered by the root rank
    integer(i32), allocatable, intent(out) :: displacements(:)

    integer(i32) :: i

    allocate( displacements, mold = receive_counts )
    displacements(1) = 0
    do i = 2, size( receive_counts )
      displacements(i) = sum( receive_counts(1 : i - 1) )
    end do
  end subroutine 

end module