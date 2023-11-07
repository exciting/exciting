module rttddft_io_serial
  use errors_warnings, only: terminate_if_false
  use m_getunit, only: getunit
  use mod_mpi_env, only: mpiinfo
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
  subroutine read_array_rank4( file_name, first, array, mpi_env )
    character(len=*), intent(in) :: file_name
    integer(i32), intent(in) :: first
    complex(dp), intent(out) :: array(:, :, :, first:)
    type(mpiinfo), intent(inout) :: mpi_env
    
    integer(i32) :: i, last, unit, iostat
    integer(long_int) :: size_block

    last = ubound( array, 4 )
    call getunit( unit )
    inquire( ioLength=size_block ) array(:, :, :, first)
    open( unit, file=trim(file_name), action='READ', &
      & form='UNFORMATTED', access='DIRECT', recl=size_block, iostat=iostat )
    call terminate_if_false( mpi_env, iostat == 0, "Error opening file: " //trim(file_name) )
    do i = first, last
      read( unit, rec=i ) array(:, :, :, i)
    end do
    close( unit )
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
    type(mpiinfo), intent(inout), optional :: mpi_env
    
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call read_array_rank4( file_name, lbound(ptr_rank4, 4), ptr_rank4, mpi_env )
  end subroutine

  subroutine write_array_rank4( file_name, first, array, mpi_env )
    character(len=*), intent(in) :: file_name
    integer(i32), intent(in) :: first
    complex(dp), intent(in) :: array(:, :, :, first:)
    type(mpiinfo), intent(inout) :: mpi_env
    
    integer(i32) :: i, last, unit, iostat
    integer(long_int) :: size_block

    last = ubound( array, 4 )
    call getunit( unit )
    inquire( ioLength=size_block ) array(:, :, :, first)
    open( unit, file=trim(file_name), action='WRITE',&
        & form='UNFORMATTED', access='DIRECT', recl=size_block, iostat=iostat)
    call terminate_if_false( mpi_env, iostat == 0, "Error opening file: " //trim(file_name) )
    do i = first, last
      write( unit, rec=i ) array(:, :, :, i)
    end do
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
    type(mpiinfo), intent(inout), optional :: mpi_env
    ! Local variables
    complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

    call map_array_to_pointer( first, array, ptr_rank4 )
    call write_array_rank4( file_name, lbound(ptr_rank4, 4), ptr_rank4, mpi_env )
  end subroutine

end module