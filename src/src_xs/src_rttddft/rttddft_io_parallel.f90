#ifdef MPI
!> Module for reading/writing binary files in parallel mode (using MPI subroutines)
module rttddft_io_parallel
  use asserts, only: assert
  use math_utils, only: all_close
  use mod_mpi_env, only: mpiinfo
  ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault with openmpi
  use mpi_f08, only: mpi_file_open, mpi_file_iread_at, mpi_file_iwrite_at, mpi_wait, mpi_f_sync_reg, &
    MPI_ASYNC_PROTECTS_NONBLOCKING, MPI_COMM, MPI_DATATYPE, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, &
    MPI_FILE, MPI_INFO_NULL, MPI_INTEGER, MPI_MODE_CREATE, MPI_MODE_RDONLY, MPI_MODE_WRONLY, &
    MPI_OFFSET_KIND, MPI_REQUEST, MPI_REQUEST_NULL, MPI_STATUS, MPI_SUCCESS
  use modmpi, only: terminate_if_false
  use precision, only: i32, dp
  use rttddft_arrays_utils, only: map_array_to_pointer
  
  implicit none
  
  private
  
  integer(i32), parameter :: bytes_complex_dp = sizeof( cmplx(0_dp, 0_dp, dp) )
  integer(i32), parameter :: bytes_real_dp = sizeof( real(0_dp, dp) )
  integer(i32), parameter :: bytes_int_i32 = sizeof( int(0, i32) )

  public :: read_array, write_array

  enum, bind(C)
    enumerator :: io_mode
    enumerator :: read_mode, write_mode
  end enum

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
    module procedure :: n_bytes_cmplx_dp
    module procedure :: n_bytes_real_dp
    module procedure :: n_bytes_int_i32
  end interface

  interface mpi_write_data
    module procedure :: mpi_write_data_complex_dp
    module procedure :: mpi_write_data_real_dp
    module procedure :: mpi_write_data_integer_i32
    module procedure :: mpi_write_data_with_header
  end interface

  interface mpi_read_data
    module procedure :: mpi_read_data_complex_dp
    module procedure :: mpi_read_data_real_dp
    module procedure :: mpi_read_data_integer_i32
    module procedure :: mpi_read_data_with_header
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

  function n_bytes_real_dp(array) result(bytes)
    !> array to measure the number of bytes
    real(dp), intent(in) :: array(..)
    integer(i32) :: bytes
    bytes = size(array) * bytes_real_dp
  end function

  function n_bytes_int_i32(array) result(bytes)
    !> array to measure the number of bytes
    integer(i32), intent(in) :: array(..)
    integer(i32) :: bytes
    bytes = size(array) * bytes_int_i32
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
    integer(MPI_OFFSET_KIND) :: bytes_block_
    bytes_block_ = bytes_block
    allocate( offset(first:last), source=[((i-1)*bytes_block_, i = first, last)] )
  end subroutine

  !> Read an array of rank=4 by chuncks
  subroutine read_array_rank4( file_name, first, array, descriptors, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 4th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to read from binary file
    complex(dp), contiguous, intent(out) :: array(:, :, :, first:)
    !> If present, these descriptors must match those stored in the file
    real(dp), contiguous, optional, intent(in) :: descriptors(:, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(in):: mpi_env    
    
    integer(i32) :: i, ierr, bytes
    integer(i32), allocatable :: dims(:, :)
    integer(MPI_OFFSET_KIND), allocatable :: offset(:)
    real(dp), allocatable :: descriptors_in_file(:, :)
    type(MPI_FILE) :: unit
    real(dp), parameter :: tol = 1.0e-7_dp

    bytes = n_bytes( array(:, :, :, first) )
    associate( a => size( array, 1 ), b => size( array, 2 ), c => size( array, 3 ), last => ubound( array, 4 ) )
      if( present( descriptors ) ) then
        call assert( ubound( descriptors, 2 ) == last, "descriptors has incompatible dim. with array")
        allocate( descriptors_in_file(size( descriptors, 1 ), first:last), dims(3, first:last) )
        bytes = bytes + n_bytes( descriptors(:, first) ) + 3*bytes_int_i32
      end if
      call file_offset( first, last, bytes, offset )
      call mpi_open_file( file_name, mpi_env, unit, read_mode )
      do i = first, last
        if( present( descriptors ) ) then                  
          call mpi_read_data( unit, offset(i), descriptors_in_file(:, i), dims(:, i), array(:, :, :, i) )
        else
          call mpi_read_data( unit, offset(i), array(:, :, :, i) )
        end if
      end do
      call mpi_file_close( unit, ierr )
      if( present( descriptors ) ) then
        call terminate_if_false( all_close( descriptors_in_file, descriptors, tol ), "Descriptors is incongruent with what is stored in file " // file_name )
        call terminate_if_false( all( dims == spread( [a, b, c], dim=2, ncopies=size( array, 4 ) ) ), "Dims is incongruent with what is stored in file " // file_name )
      end if
    end associate
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
    type(mpiinfo), intent(in):: mpi_env    
    ! Local variables
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call read_array_rank4( file_name, lbound( ptr_rank4, 4 ), ptr_rank4, mpi_env=mpi_env )
  end subroutine

  !> Write an array of rank=4 by chuncks, including headers if `descriptors` is present
  subroutine write_array_rank4( file_name, first, array, descriptors, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 4th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to be written to binary file
    complex(dp), contiguous, intent(inout) :: array(:, :, :, first:)
    !> descriptors to be written as part of the header
    real(dp), contiguous, optional, intent(inout) :: descriptors(:, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(in):: mpi_env 

    integer(i32) :: i, bytes, ierr, dims(3)
    integer(MPI_OFFSET_KIND), allocatable :: offset(:)
    type(MPI_FILE) :: unit

    bytes = n_bytes( array(:, :, :, first) )
    associate( last => ubound( array, 4 ) )
      if( present( descriptors ) ) then
        call assert( ubound( descriptors, 2 ) == last, "descriptors has incompatible dim. with array")
        bytes = bytes + n_bytes( descriptors(:, first) ) + 3*bytes_int_i32
        dims = [size(array, 1), size(array, 2), size(array, 3)]
      end if
      call file_offset( first, last, bytes, offset )
      call mpi_open_file( file_name, mpi_env, unit, write_mode )
      do i = first, last
        if( present( descriptors ) ) then                  
          call mpi_write_data( unit, offset(i), descriptors(:, i), dims, array(:, :, :, i) )
        else
          call mpi_write_data( unit, offset(i), array(:, :, :, i) )
        end if
      end do
      call mpi_file_close( unit, ierr )
    end associate
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
    type(mpiinfo), intent(in):: mpi_env 
    ! Local variables
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)
    
    call map_array_to_pointer( first, array, ptr_rank4 )
    call write_array_rank4( file_name, lbound( ptr_rank4, 4 ), ptr_rank4, mpi_env=mpi_env )
  end subroutine

  ! MPI-IO wrappers (private)
  subroutine mpi_read_data_complex_dp( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    complex(dp), contiguous, asynchronous, intent(out) :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status
    
    request = MPI_REQUEST_NULL
    call mpi_file_iread_at( unit, offset, data_block, size(data_block), MPI_DOUBLE_COMPLEX, request, ierr )
    call mpi_wait( request, status, ierr )
    ! Ensure the correct treatment of buffers passed to nonblocking MPI-routines
    if( .not. MPI_ASYNC_PROTECTS_NONBLOCKING ) call mpi_f_sync_reg( data_block )
  end subroutine

  subroutine mpi_read_data_real_dp( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    real(dp), contiguous, asynchronous, intent(out) :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status
    
    request = MPI_REQUEST_NULL
    call mpi_file_iread_at( unit, offset, data_block, size(data_block), MPI_DOUBLE_PRECISION, request, ierr )
    call mpi_wait( request, status, ierr )
    if( .not. MPI_ASYNC_PROTECTS_NONBLOCKING ) call mpi_f_sync_reg( data_block )
  end subroutine

  subroutine mpi_read_data_integer_i32( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    integer(i32), contiguous, asynchronous, intent(out) :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status
    
    request = MPI_REQUEST_NULL
    call mpi_file_iread_at( unit, offset, data_block, size(data_block), MPI_INTEGER, request, ierr )
    call mpi_wait( request, status, ierr )
    if( .not. MPI_ASYNC_PROTECTS_NONBLOCKING ) call mpi_f_sync_reg( data_block )
  end subroutine

  !> Read a series of data: descriptor_block, dims_block, data_block
  subroutine mpi_read_data_with_header( unit, offset, descriptor_block, dims_block, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    real(dp), contiguous, asynchronous, intent(out) :: descriptor_block(..)
    integer(i32), contiguous, asynchronous, intent(out) :: dims_block(..)
    complex(dp), contiguous, asynchronous, intent(out)  :: data_block(..)

    call mpi_read_data( unit, offset, descriptor_block )
    call mpi_read_data( unit, offset + n_bytes(descriptor_block), dims_block )
    call mpi_read_data( unit, offset + n_bytes(descriptor_block) + n_bytes(dims_block), data_block )
  end subroutine

  subroutine mpi_write_data_complex_dp( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    !> N.B. intent(inout) is required by `call mpi_f_sync_reg( data_block )`
    complex(dp), contiguous, asynchronous, intent(inout)  :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status

    request = MPI_REQUEST_NULL
    call mpi_file_iwrite_at( unit, offset, data_block, size(data_block), MPI_DOUBLE_COMPLEX, request, ierr )
    call mpi_wait( request, status, ierr )
    ! Ensure the correct treatment of buffers passed to nonblocking MPI-routines
    if( .not. MPI_ASYNC_PROTECTS_NONBLOCKING ) call mpi_f_sync_reg( data_block )
  end subroutine

  subroutine mpi_write_data_real_dp( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    real(dp), contiguous, asynchronous, intent(inout)  :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status

    request = MPI_REQUEST_NULL
    call mpi_file_iwrite_at( unit, offset, data_block, size(data_block), MPI_DOUBLE_PRECISION, request, ierr )
    call mpi_wait( request, status, ierr )
    ! Ensure the correct treatment of buffers passed to nonblocking MPI-routines
    if( .not. MPI_ASYNC_PROTECTS_NONBLOCKING ) call mpi_f_sync_reg( data_block )
  end subroutine

  subroutine mpi_write_data_integer_i32( unit, offset, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    integer(i32), contiguous, asynchronous, intent(inout)  :: data_block(..)

    integer(i32) :: ierr
    type(MPI_REQUEST) :: request
    type(MPI_STATUS) :: status

    request = MPI_REQUEST_NULL
    call mpi_file_iwrite_at( unit, offset, data_block, size(data_block), MPI_INTEGER, request, ierr )
    call mpi_wait( request, status, ierr )
    ! Ensure the correct treatment of buffers passed to nonblocking MPI-routines
    if( .not. MPI_ASYNC_PROTECTS_NONBLOCKING ) call mpi_f_sync_reg( data_block )
  end subroutine

  !> Write a series of data: descriptor_block, dims_block, data_block
  subroutine mpi_write_data_with_header( unit, offset, descriptor_block, dims_block, data_block )
    type(MPI_FILE), intent(in) :: unit
    integer(MPI_OFFSET_KIND) :: offset
    real(dp), contiguous, asynchronous, intent(inout) :: descriptor_block(..)
    integer(i32), contiguous, asynchronous, intent(inout) :: dims_block(..)
    complex(dp), contiguous, asynchronous, intent(inout)  :: data_block(..)

    call mpi_write_data( unit, offset, descriptor_block )
    call mpi_write_data( unit, offset + n_bytes(descriptor_block), dims_block )
    call mpi_write_data( unit, offset + n_bytes(descriptor_block) + n_bytes(dims_block), data_block )
  end subroutine

  subroutine mpi_open_file( file_name, mpi_env, unit, mode )
    character(len=*), intent(in) :: file_name
    type(mpiinfo), intent(in) :: mpi_env
    type(MPI_FILE), intent(out) :: unit
    integer(kind(io_mode)), intent(in) :: mode

    integer(i32) :: ierr
    type(MPI_COMM) :: handle

    call assert( mode==read_mode .or. mode==write_mode, 'mode must be read or write' )

    handle%mpi_val = mpi_env%comm
    select case( mode ) 
      case( read_mode )
        call mpi_file_open( handle, trim(file_name), &
          MPI_MODE_RDONLY, MPI_INFO_NULL, unit, ierr )
      case( write_mode )
        call mpi_file_open( handle, trim(file_name), &
          MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, unit, ierr )
    end select
    call terminate_if_false( ierr == MPI_SUCCESS, "Error opening file "//trim(file_name) )
  end subroutine

end module 
#endif
