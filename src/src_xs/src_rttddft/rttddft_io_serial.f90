! MRM (2025) For Cray compiler we need 64 bits integers. This have 
! the side effect on affecting the logical default to 64 bit.
! Therefore, we need to indicate the size of logical.
module rttddft_io_serial
  use asserts, only: assert
  use math_utils, only: all_close
  use mod_mpi_env, only: mpiinfo
  use modmpi, only: terminate_if_false
  use precision, only: i32, long_int, dp
  use rttddft_arrays_utils, only: map_array_to_pointer

  implicit none
  
  private
  
  public :: read_array, write_array

  interface read_array
    module procedure :: read_array_rank4
    module procedure :: read_array_rank5
  end interface

  interface write_array
    module procedure :: write_array_rank4
    module procedure :: write_array_rank5
  end interface

contains
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
    
    integer(i32) :: i, unit, iostat
    integer(i32), allocatable :: dims(:, :)
    integer(long_int) :: size_block
    real(dp), allocatable :: descriptors_in_file(:, :)
    real(dp), parameter :: tol = 1.0e-7_dp

    associate( m => size(array, 1, kind=i32), n => size(array, 2, kind=i32), &
               k => size(array, 3, kind=i32), last => ubound(array, 4, kind=i32) )
      if( present(descriptors) ) then
        call assert( logical(ubound( descriptors, 2, kind=i32 ) == last, kind=i32) , "incompatible descriptors and array")
        allocate( descriptors_in_file(size(descriptors, 1), first:last), dims(3, first:last) )
        inquire( ioLength=size_block ) descriptors(:, first), dims(:, first), array(:, :, :, first)
      else
        inquire( ioLength=size_block ) array(:, :, :, first)
      end if
      open( newunit=unit, file=trim(file_name), action='READ', form='UNFORMATTED', access='DIRECT', recl=size_block, iostat=iostat )
      call terminate_if_false( logical(iostat == 0, kind=i32), "Error opening file: " //trim(file_name) )
      do i = first, last
        if( present(descriptors) ) then
          read( unit, rec=i ) descriptors_in_file(:, i), dims(:, i), array(:, :, :, i)               
        else
          read( unit, rec=i ) array(:, :, :, i)
        end if
      end do
      close( unit )
      if( present(descriptors) ) then
        call terminate_if_false( all_close(descriptors_in_file, descriptors, tol), "Descriptors is incongruent with what is stored in file " // file_name )
        call terminate_if_false( all( logical(dims == spread([m, n, k], dim=2, ncopies=size(array, 4, kind=i32)), kind=i32)), &
                                 "Dims is incongruent with what is stored in file " // file_name )
      end if
    end associate
  end subroutine

  !> Remap an array of rank 5 to an array of rank 4 using a pointer
  !> and then use [[read_array_rank4]]
  subroutine read_array_rank5( file_name, first, array, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 5th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to be read from binary file
    complex(dp), contiguous, target, intent(out) :: array(:, :, :, :, first:)
    !> MPI environment. Ignored, but needed to have the same interface as its 
    !> parallel counterpart in [[rttddft_io_parallel]]
    type(mpiinfo), intent(in), optional :: mpi_env
    
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call read_array_rank4( file_name, lbound(ptr_rank4, 4, kind=i32), ptr_rank4, mpi_env=mpi_env )
  end subroutine

  !> Write an array of rank=4 by chuncks, including headers if `descriptors` is present
  subroutine write_array_rank4( file_name, first, array, descriptors, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 4th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to be written to binary file
    complex(dp), contiguous, intent(in) :: array(:, :, :, first:)
    !> descriptors to be written before each chunck
    real(dp), contiguous, optional, intent(in) :: descriptors(:, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(in):: mpi_env 

    integer(i32) :: i, unit, iostat
    integer(long_int) :: size_block

    associate( m => size(array, 1, kind=i32), n => size(array, 2, kind=i32), &
               k => size(array, 3, kind=i32), last => ubound( array, 4, kind=i32) )
      if( present(descriptors) ) then
        call assert( logical(ubound( descriptors, 2, kind=i32 ) == last, kind=i32), "incompatible descriptors and array")
        inquire( ioLength=size_block ) descriptors(:, first), m, n, k, array(:, :, :, first)
      else
        inquire( ioLength=size_block ) array(:, :, :, first)
      end if
      open( newunit=unit, file=trim(file_name), action='WRITE', form='UNFORMATTED', access='DIRECT', recl=size_block, iostat=iostat)
      call terminate_if_false( logical(iostat == 0, kind=i32), "Error opening file: " //trim(file_name) )
      do i = first, last
        if( present(descriptors) ) then
          write( unit, rec=i ) descriptors(:, i), m, n, k, array(:, :, :, i)
        else
          write( unit, rec=i ) array(:, :, :, i)
        end if
      end do
    end associate
    close( unit )
  end subroutine

  !> Remap an array of rank 5 to an array of rank 4 using a pointer
  !> and then use [[write_array_rank4]]
  subroutine write_array_rank5( file_name, first, array, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 5th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> Array to be written to binary file
    complex(dp), contiguous, target, intent(in) :: array(:, :, :, :, first:)
    !> MPI environment. Ignored, but needed to have the same interface as its 
    !> parallel counterpart in [[rttddft_io_parallel]]
    type(mpiinfo), intent(in), optional :: mpi_env
    ! Local variables
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call write_array_rank4( file_name, lbound(ptr_rank4, 4, kind=i32), ptr_rank4, mpi_env=mpi_env )
  end subroutine

end module
