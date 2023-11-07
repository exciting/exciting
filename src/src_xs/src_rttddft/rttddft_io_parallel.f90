#ifdef MPI
!> Module for reading/writing binary files in parallel mode (using MPI subroutines)
module rttddft_io_parallel
  use asserts, only: assert
  use errors_warnings, only: terminate_if_false
  use mod_mpi_env, only: mpiinfo
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault with openmpi
  use mpi_f08, only: mpi_file_open, mpi_file_iread_at, mpi_file_iwrite_at, mpi_wait, &
    MPI_DATATYPE, MPI_COMM, MPI_DOUBLE_COMPLEX, MPI_FILE, MPI_REQUEST, MPI_OFFSET_KIND, &
    MPI_STATUS, MPI_INFO_NULL, MPI_MODE_CREATE, MPI_MODE_RDONLY, MPI_MODE_WRONLY, &
    MPI_SUCCESS
  use precision, only: i32, dp
  use rttddft_arrays_utils, only: map_array_to_pointer
  
  implicit none
  
  private
  
  integer(i32), parameter :: bytes_complex_dp = sizeof( cmplx(0_dp,0_dp,dp) )

  public :: read_array, write_array

  ! This can be expanded to more ranks when needed
  interface read_array
    module procedure :: read_array_rank4
    module procedure :: read_array_rank5
  end interface

  ! This can be expanded to more ranks when needed
  interface write_array
    module procedure :: write_array_rank4
    module procedure :: write_array_rank5
  end interface

  interface n_bytes
    module procedure n_bytes_cmplx_dp
  end interface

contains
  !> Return the number of bytes of a complex array
  !> This could be overloaded w.r.t. type
  function n_bytes_cmplx_dp(array) result(bytes)
    !> array to measure the number of bytes
    complex(dp), intent(in) :: array(..)
    integer(i32) :: bytes
    bytes = size(array) * bytes_complex_dp
  end function
  
  !> Return an array of offsets. Convention: the first chunck (with `first=0`) 
  !> should be placed at position 0 (`offset=0`)
  subroutine file_offset(first, last, bytes_block, offset)
    !> first index of array that will be written/read
    integer(i32), intent(in) :: first
    !> last index of array
    integer(i32), intent(in) :: last
    !> number of bytes of each block (chunck of an array)
    integer(i32), intent(in) :: bytes_block 
    !> array with the offsets
    integer(MPI_OFFSET_KIND), allocatable, intent(out) :: offset(:)

    integer(i32) :: i
    allocate( offset(first:last), source=int([((i-1)*bytes_block, i = first, last)], MPI_OFFSET_KIND) )
  end subroutine

  !> Read an array of rank=4 by chuncks
  subroutine read_array_rank4( file_name, first, array, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 4th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to be read from binary file
    complex(dp), intent(out) :: array(:, :, :, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(inout):: mpi_env    
    
    integer(i32) :: i, last, ierr
    integer(MPI_OFFSET_KIND), allocatable :: offset(:)
    type(MPI_FILE) :: unit
    
    last = ubound( array, 4 )
    call file_offset( first, last, n_bytes( array(:,:,:,first) ), offset )
    call mpi_open_file( file_name, mpi_env, unit, 'read' )
    do i = first, last                        
      call mpi_read_data( unit, offset(i), array(:, :, :, i) )
    end do
    call mpi_file_close( unit, ierr )
  end subroutine

  !> Remap an array of rank 5 to an array of rank 4 using a pointer
  subroutine read_array_rank5( file_name, first, array, mpi_env )
    !> Name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> First index along 5th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> Array to be read from binary file
    complex(dp), contiguous, target, intent(out) :: array(:, :, :, :, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(inout):: mpi_env    
    ! Local variables
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call read_array_rank4( file_name, lbound( ptr_rank4, 4 ), ptr_rank4, mpi_env )
  end subroutine

  !> Write an array of rank=4 by chuncks
  subroutine write_array_rank4( file_name, first, array, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 4th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to be written to binary file
    complex(dp), intent(in) :: array(:, :, :, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(inout):: mpi_env 

    integer(i32) :: i, last, ierr
    integer(MPI_OFFSET_KIND), allocatable :: offset(:)
    type(MPI_FILE) :: unit
    
    last = ubound( array, 4 )
    call file_offset( first, last, n_bytes( array(:,:,:,first) ), offset )
    call mpi_open_file( file_name, mpi_env, unit, 'write' )
    do i = first, last    
      call mpi_write_data( unit, offset(i), array(:, :, :, i) )
    end do    
    call mpi_file_close( unit, ierr )
  end subroutine

  !> Write an array of rank=5 by chuncks
  subroutine write_array_rank5( file_name, first, array, mpi_env )
    !> Name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> First index along 5th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> Array to be written to binary file
    complex(dp), contiguous, target, intent(in) :: array(:, :, :, :, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(inout):: mpi_env 
    ! Local variables
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)
    
    call map_array_to_pointer( first, array, ptr_rank4 )
    call write_array_rank4( file_name, lbound( ptr_rank4, 4 ), ptr_rank4, mpi_env )
  end subroutine

  ! MPI-IO wrappers (they are private)
  subroutine mpi_read_data( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    complex(dp), intent(out) :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status
    call mpi_file_iread_at( unit, offset, data_block, &
        size(data_block), MPI_DOUBLE_COMPLEX, request, ierr )
    call mpi_wait( request, status, ierr )
  end subroutine

  subroutine mpi_write_data( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    complex(dp), intent(in)  :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status
    call mpi_file_iwrite_at( unit, offset, data_block, &
        size(data_block), MPI_DOUBLE_COMPLEX, request, ierr )
    call mpi_wait( request, status, ierr )
  end subroutine

  subroutine mpi_open_file( file_name, mpi_env, unit, mode )
    character(len=*), intent(in) :: file_name
    type(mpiinfo), intent(inout) :: mpi_env
    type(MPI_FILE), intent(out) :: unit
    character(len=*), intent(in) :: mode

    integer(i32) :: ierr
    type(MPI_COMM) :: handle

    call assert( trim(mode)=='read' .or. trim(mode)=='write', 'mode must be read or write' )

    handle%mpi_val = mpi_env%comm
    select case( trim(mode) ) 
      case('read')
        call mpi_file_open( handle, trim(file_name), &
          MPI_MODE_RDONLY, MPI_INFO_NULL, unit, ierr )
      case('write')
        call mpi_file_open( handle, trim(file_name), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, unit, ierr )
    end select
    call terminate_if_false( mpi_env, ierr == MPI_SUCCESS, "Error opening file"//trim(file_name) )
  end subroutine

end module 
#endif