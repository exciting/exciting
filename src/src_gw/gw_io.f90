!> This module is designed to centralize IO operation in GW
module gw_io
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use m_getunit, only: getunit
  use precision, only: dp, i32

  implicit none

  private

  !> Unit used for the general GW output file (`GW_INFO.OUT`)
  integer(i32), public, protected :: fgw
  !> Default name of the general GW output file
  character(len=*), parameter :: filename_gwinfo = 'GW_INFO.OUT'
  !> Default extension
  character(len=*), parameter :: default_file_extension = '.OUT'
  !> Accepted file format = 'text'
  character(len=*), parameter :: file_format_text = 'text'
  !> Accepted file format = 'binary'
  character(len=*), parameter :: file_format_binary = 'binary'
  
  public :: open_gwinfo, write_to_gwinfo, write_to_gwinfo_boxmessage, &
    build_file_name, write_to_file, read_from_file
  
  interface write_to_file
    module procedure write_matrix_to_file
    module procedure write_matrix_to_file_given_lbounds
    module procedure write_tensor_of_rank_3_to_file
    module procedure write_tensor_of_rank_3_to_file_given_lbounds
  end interface
  
  interface build_file_name
    module procedure build_file_name_with_integer
    module procedure build_file_name_only_adding_extension
  end interface

  interface read_from_file
    module procedure read_matrix_from_file
    module procedure read_tensor_of_rank_3_from_file
  end interface

contains

!> Open the `GW_INFO.OUT`
subroutine open_gwinfo( )
  
  call open_text_file( filename_gwinfo , 'write', fgw )

end subroutine


!> Write a string into `GW_INFO.OUT`
subroutine write_to_gwinfo( string )
  character(len=*), intent(in) :: string

  write( fgw, * ) string

end subroutine


!> Write a string surrounded by a box of characters into `GW_INFO.OUT`
subroutine write_to_gwinfo_boxmessage( char, string )
  character, intent(in) :: char
  character(len=*), intent(in) :: string

  call BoxMSG( fgw, char, string )

end subroutine


!> Append an integer to a string together with the default extension, and 
!> store the result in `file_name`
subroutine build_file_name_with_integer( base_name, int, file_name )
  !> String containing the base name, to which an integer will be added
  character(len=*), intent(in)  :: base_name
  !> Integer to be appended to the string `base_name`
  integer(i32), intent(in)      :: int
  !> Resulting string from appending `int` to the the end of `base_name`, and adding the extension
  character(len=*), intent(out) :: file_name

  write( file_name, * ) int
  file_name = trim( adjustl(base_name) ) // trim( adjustl( file_name ) ) // default_file_extension

end subroutine


!> Add the default extension to `base_name` and store the result in `file_name`
subroutine build_file_name_only_adding_extension( base_name, file_name )
  !> String containing the base name, to which the default extension will be added
  character(len=*), intent(in)  :: base_name
  !> Resulting string from adding the default extension to `base_name`
  character(len=*), intent(out) :: file_name

  file_name = trim( adjustl(base_name) ) // default_file_extension

end subroutine


!> Write a matrix (array of rank 2) to a file
subroutine write_matrix_to_file( file_name, matrix, file_format )
  character(len=*), intent(in) :: file_name
  complex(dp), intent(in) :: matrix(:, :)
  character(len=*), intent(in) :: file_format

  call write_matrix_to_file_given_lbounds( file_name, matrix, [1, 1], file_format )

end subroutine


!> Write a matrix (array of rank 2) to a file, given its lower bounds
subroutine write_matrix_to_file_given_lbounds( file_name, matrix, lbounds, file_format )
  !> File name where to write
  character(len=*), intent(in) :: file_name
  !> lbounds of `matrix`
  integer(i32), intent(in) :: lbounds(2)
  !> Matrix to be written into the file
  complex(dp), intent(in) :: matrix(lbounds(1):, lbounds(2):)
  !> File format of output
  character(len=*), intent(in) :: file_format
  
  integer(i32) :: unit, i

  call open_file_generic( file_name, 'write', file_format, unit )
  call write_header_to_file( unit, matrix, lbounds )
  do i = lbounds(2), ubound( matrix, 2 )
    write( unit, * ) matrix(:, i)
  end do
  close( unit )
end subroutine
  

!> Write a tensor of rank 3 to a file
subroutine write_tensor_of_rank_3_to_file( file_name, tensor, file_format )
  !> File name where to write
  character(len=*), intent(in) :: file_name
  !> Tensor to be written into the file
  complex(dp), intent(in) :: tensor(:, :, :)
  !> Format of the output file
  character(len=*), intent(in) :: file_format

  call write_tensor_of_rank_3_to_file_given_lbounds( file_name, tensor, [1, 1, 1], file_format )

end subroutine


subroutine write_tensor_of_rank_3_to_file_given_lbounds( file_name, tensor, lbounds, file_format )
  character(len=*), intent(in) :: file_name
  integer(i32), intent(in) :: lbounds(3)
  complex(dp), intent(in) :: tensor(lbounds(1):, lbounds(2):, lbounds(3):)
  character(len=*), intent(in) :: file_format

  integer(i32) :: unit, i, j

  call open_file_generic( file_name, 'write', file_format, unit )
  call write_header_to_file( unit, tensor, lbound(tensor) )
  do i = lbound( tensor, 3 ), ubound( tensor, 3 )
    do j = lbound( tensor, 2 ), ubound( tensor, 2 )
      write( unit, * ) tensor(:, j, i)
    end do
  end do
  close( unit )

end subroutine


!> (private) Write a header to an output file
!> The header contains the rank of the array, its lbounds and ubounds
subroutine write_header_to_file( unit, array, lbounds )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Array to be written into the file
  complex(dp), intent(in) :: array(..)
  !> Lbounds of `array`
  integer(i32), intent(in) :: lbounds(:)

  call assert( size(lbounds) == rank(array), 'Incompatible array rank and lbounds' )
  write( unit, * ) rank( array )
  write( unit, * ) lbounds, ubound( array ) + lbounds - 1

end subroutine


!> (private) Read the header of an input file
!> The header contains the rank of the array, its lbounds and ubounds
subroutine read_header_of_file( unit, rank_of_array, lbounds, ubounds )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Rank of the array stored in the file
  integer(i32), intent(out) :: rank_of_array
  !> Lbounds of the array
  integer(i32), intent(out) :: lbounds(:)
  !> Ubounds of the array
  integer(i32), intent(out) :: ubounds(:)
  
  call assert( size(lbounds) == size(ubounds), 'lbounds and ubounds must have same size')
  read( unit, * ) rank_of_array
  call terminate_if_false( rank_of_array==size(lbounds), 'Incompatible rank of array' )
  read( unit, * ) lbounds, ubounds

end subroutine


!> Check if a file exists and terminates the execution if it does not
subroutine terminate_if_file_does_not_exist( file_name )
  !> Name of the file to check
  character(len=*), intent(in) :: file_name

  logical :: ok
  inquire( file=trim(file_name), exist=ok )
  call terminate_if_false( ok, 'File '//trim(file_name)//' not found')

end subroutine
  
  
!> Read a matrix (array of rank 2) from a file
subroutine read_matrix_from_file( file_name, matrix, file_format )
  !> Name of the file to read
  character(len=*), intent(in)          :: file_name
  !> Matrix to store the data read from the file
  complex(dp), intent(out), allocatable :: matrix(:, :)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: unit, i, read_rank, lbounds(2), ubounds(2)
  integer(i32), parameter :: expected_rank = 2 !rank of a matrix

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, read_rank, lbounds, ubounds )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no matrix' )
  allocate( matrix(lbounds(1):ubounds(1), lbounds(2):ubounds(2)) )
  do i = lbounds(2), ubounds(2)
    read( unit, * ) matrix(:, i)
  end do
  close( unit )

end subroutine
  

!> Read a tensor of rank 3 stored in a file
subroutine read_tensor_of_rank_3_from_file( file_name, tensor, file_format )
  !> Name of the file to read
  character(len=*), intent(in)          :: file_name
  !> Tensor to store the data read from the file
  complex(dp), intent(out), allocatable :: tensor(:, :, :)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: unit, i, j, read_rank, lbounds(3), ubounds(3)
  integer(i32), parameter :: expected_rank = 3 !rank of a tensor

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, read_rank, lbounds, ubounds )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no tensor of rank 3' )
  allocate( tensor(lbounds(1):ubounds(1), lbounds(2):ubounds(2), lbounds(3):ubounds(3)) )
  do i = lbounds(3), ubounds(3) 
    do j = lbounds(2), ubounds(2) 
      read( unit, * ) tensor(:, j, i)
    end do
  end do
  close( unit )

end subroutine


!> (private) generic subroutine to open a file
subroutine open_file_generic( file_name, action, file_format, unit )
  !> name of file to open
  character(len=*), intent(in) :: file_name
  !> action (typically: read/write)
  character(len=*), intent(in) :: action
  !> format: currently, only binary and text are accepted
  character(len=*), intent(in) :: file_format
  !> unit number of file to open
  integer, intent(out) :: unit

  call assert( trim(file_format)==file_format_text .or. &
    trim(file_format)==file_format_binary, &
    'file_format must be '// file_format_text // ' or ' // file_format_binary )
  
  select case( trim( file_format ) )
    case( file_format_text )
      call open_text_file( file_name, action, unit )
    case( file_format_binary )
      call open_binary_file( file_name, action, unit )
  end select
end subroutine


!> Open a text file
subroutine open_text_file( file_name, action, unit )
  !> Name of the file
  character(len=*), intent(in)  :: file_name
  !> Action: should be e. g. "read", "write"
  character(len=*), intent(in)  :: action
  !> Unit of the file (that will be opened)
  integer(i32), intent(out)     :: unit

  call getunit( unit )
  open( unit, file=trim(file_name), action=action, form='formatted' )

end subroutine

  
!> Open a binary file
subroutine open_binary_file( file_name, action, unit )
  !> Name of the file
  character(len=*), intent(in)  :: file_name
  !> Action: should be e. g. "read", "write"
  character(len=*), intent(in)  :: action
  !> Unit of the file (that will be opened)
  integer(i32), intent(out)     :: unit

  call getunit( unit )
  open( unit, file=trim(file_name), action=action, form='unformatted', access='stream' )
  
end subroutine
  
end module
