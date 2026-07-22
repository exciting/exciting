module rttddft_io_serial
#include "asserts.fpp"
  use math_utils, only: all_close
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use mod_mpi_env, only: mpiinfo
  use modmpi, only: terminate_if_false
  use precision, only: i32, long_int, dp
  use rttddft_arrays_utils, only: map_array_to_pointer

  implicit none
  
  private
  
  public :: read_array, read_three_arrays, write_array, write_three_arrays

  interface read_array
    module procedure :: read_array_rank3
    module procedure :: read_array_rank4
    module procedure :: read_array_rank5
  end interface

  interface read_three_arrays
    module procedure :: read_three_arrays_rank3
  end interface

  interface write_array
    module procedure :: write_array_rank4
    module procedure :: write_array_rank5
  end interface

  interface write_three_arrays
    module procedure :: write_three_arrays_rank3
  end interface

contains
  !> Read an array of rank=3 by chuncks
  subroutine read_array_rank3( file_name, first, array, descriptors, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 4th dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> array to read from binary file
    complex(dp), contiguous, intent(out) :: array(:, :, first:)
    !> If present, these descriptors must match those stored in the file
    real(dp), contiguous, optional, intent(in) :: descriptors(:, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(in):: mpi_env    
    
    integer(i32) :: i, unit
    integer(i32), allocatable :: dims(:, :)
    integer(long_int) :: size_block
    real(dp), allocatable :: descriptors_in_file(:, :)
    real(dp), parameter :: tol = 1.0e-7_dp

    associate( m => size(array, 1), n => size(array, 2), &
               last => ubound(array, 3) )
      if( present(descriptors) ) then
        CALL_ASSERT( logical(ubound( descriptors, 2 ) == last) , "incompatible descriptors and array")
        allocate( descriptors_in_file(size(descriptors, 1), first:last), dims(1, first:last) )
        ! Trick: dims is a 1xN array to use an existing interface of inquire_large
        call inquire_large( size_block, descriptors(:, first), dims(:, first), array(:, :, first) )
      else
        call inquire_large( size_block, array(:, :, first) )
      end if
      call open_direct_unformatted_large( unit, trim(file_name), "read", size_block, "old" )
      do i = first, last
        if( present(descriptors) ) then
          read( unit, rec=i ) descriptors_in_file(:, i), dims(1, i), array(:, :, i)               
        else
          read( unit, rec=i ) array(:, :, i)
        end if
      end do
      close( unit )
      if( present(descriptors) ) then
        call terminate_if_false( all_close(descriptors_in_file, descriptors, tol), "Descriptors is incongruent with what is stored in file " // file_name )
        call terminate_if_false( all( logical(dims(1, first:last) == spread(m, dim=1, ncopies=size(array, 3)))), &
                                 "Dims is incongruent with what is stored in file " // file_name )
      end if
    end associate
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
    
    integer(i32) :: i, unit
    integer(i32), allocatable :: dims(:, :)
    integer(long_int) :: size_block
    real(dp), allocatable :: descriptors_in_file(:, :)
    real(dp), parameter :: tol = 1.0e-7_dp

    associate( m => size(array, 1), n => size(array, 2), &
               k => size(array, 3), last => ubound(array, 4) )
      if( present(descriptors) ) then
        CALL_ASSERT( logical(ubound( descriptors, 2 ) == last) , "incompatible descriptors and array")
        allocate( descriptors_in_file(size(descriptors, 1), first:last), dims(3, first:last) )
        call inquire_large( size_block, descriptors(:, first), dims(:, first), array(:, :, :, first) )
      else
        call inquire_large( size_block, array(:, :, :, first) )
      end if
      call open_direct_unformatted_large( unit, trim(file_name), "read", size_block, "old" )
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
        call terminate_if_false( all( logical(dims == spread([m, n, k], dim=2, ncopies=size(array, 4)))), &
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
    type(mpiinfo), intent(in) :: mpi_env
    
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call read_array_rank4( file_name, lbound(ptr_rank4, 4), ptr_rank4, mpi_env=mpi_env )
  end subroutine

  !> Read three arrays of rank=3 by chuncks
  subroutine read_three_arrays_rank3( file_name, first, array_1, array_2, array_3, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 3rd dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> 1st array to be read from the binary file
    complex(dp), contiguous, intent(out) :: array_1(:, :, first:)
    !> 2nd array to be read from the binary file
    complex(dp), contiguous, intent(out) :: array_2(:, :, first:)
    !> 3rd array to be read from the binary file
    complex(dp), contiguous, intent(out) :: array_3(:, :, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(in):: mpi_env 

    integer(i32) :: i, last, unit
    integer(long_int) :: size_block

    last = ubound( array_1, 3 )
    CALL_ASSERT( all( shape(array_1) == shape(array_2) ) .and. all( shape(array_1) == shape(array_3) ), "shape mismatch" )
    call inquire_large( size_block, array_1(:, :, first),array_2(:, :, first), array_3(:, :, first) )
    call open_direct_unformatted_large( unit, trim(file_name), "read", size_block, "old" )
    do i = first, last
      read( unit, rec=i ) array_1(:, :, i), array_2(:, :, i), array_3(:, :, i)
    end do
    close( unit )
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

    integer(i32) :: i, unit
    integer(long_int) :: size_block

    associate( m => size(array, 1), n => size(array, 2), &
               k => size(array, 3), last => ubound( array, 4) )
      if( present(descriptors) ) then
        CALL_ASSERT( logical(ubound( descriptors, 2 ) == last), "incompatible descriptors and array")
        call inquire_large( size_block, descriptors(:, first), [m, n, k], array(:, :, :, first) )
      else
        call inquire_large( size_block, array(:, :, :, first) )
      end if
      call open_direct_unformatted_large( unit, trim(file_name), "write", size_block, "unknown" )
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
    call write_array_rank4( file_name, lbound(ptr_rank4, 4), ptr_rank4, mpi_env=mpi_env )
  end subroutine

  !> Write three arrays of rank=3 by chuncks
  subroutine write_three_arrays_rank3( file_name, first, array_1, array_2, array_3, mpi_env )
    !> name of the file where the array is stored
    character(len=*), intent(in) :: file_name
    !> first index along 3rd dim (needed to determine offsets)
    integer(i32), intent(in) :: first
    !> 1st array to be written to the binary file
    complex(dp), contiguous, intent(in) :: array_1(:, :, first:)
    !> 2nd array to be written to the binary file
    complex(dp), contiguous, intent(in) :: array_2(:, :, first:)
    !> 3rd array to be written to the binary file
    complex(dp), contiguous, intent(in) :: array_3(:, :, first:)
    !> MPI environment. The corresponding MPI processes will read from file
    type(mpiinfo), intent(in):: mpi_env 

    integer(i32) :: i, last, unit
    integer(long_int) :: size_block

    last = ubound( array_1, 3 )
    CALL_ASSERT( all( shape(array_1) == shape(array_2) ) .and. all( shape(array_1) == shape(array_3) ), "shape mismatch" )
    call inquire_large( size_block, array_1(:, :, first),array_2(:, :, first), array_3(:, :, first) )
    call open_direct_unformatted_large( unit, trim(file_name), "write", size_block, "unknown" )
    do i = first, last
      write( unit, rec=i ) array_1(:, :, i), array_2(:, :, i), array_3(:, :, i)
    end do
    close( unit )
  end subroutine
end module
