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
  
  public :: open_gwinfo, write_to_gwinfo, write_to_gwinfo_boxmessage, &
    build_file_name, write_to_file, read_from_file
  
  interface write_to_file
    module procedure write_matrix_to_file_with_header
    module procedure write_tensor_of_rank_3_to_file
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
subroutine write_matrix_to_file_with_header( file_name, matrix, lbounds, binary_format )
  !> File name where to write
  character(len=*), intent(in) :: file_name
  !> lbounds of `matrix`
  integer(i32), intent(in) :: lbounds(2)
  !> Matrix to be written into the file
  complex(dp), intent(in) :: matrix(lbounds(1):, lbounds(2):)
  !> If true, the output has binary format
  logical, intent(in) :: binary_format

  integer(i32) :: unit, i, m, n, mini, mend, nini, nend
  integer(i32), parameter :: rank_of_array = 2 ! A matrix is an array of rank = 2

  if( binary_format ) then
    call open_binary_file( file_name, 'write', unit )
  else
    call open_text_file( file_name, 'write', unit )
  end if
  write( unit, * ) rank_of_array
  m = size( matrix, 1 )
  n = size( matrix, 2 )
  mini = lbounds(1)
  nini = lbounds(2)
  mend = mini + m -1
  nend = nini + n -1
  write( unit, * ) mini, nini, mend, nend 
  do i = nini, nend 
    write( unit, * ) matrix(mini:mend, i)
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

  integer(i32) :: unit, i, j
  integer(i32), parameter :: rank_of_array = 3

  call assert( trim(file_format)=='text' .or. trim(file_format)=='binary', &
    'file_format must be text or binary' )
  
  select case( trim( file_format ) )
    case( 'text' )
      call open_text_file( file_name, 'write', unit )
    case( 'binary' )
      call open_binary_file( file_name, 'write', unit )
  end select
  write( unit, * ) rank_of_array
  write( unit, * ) lbound( tensor ), ubound( tensor )
  do i = lbound( tensor, 3 ), ubound( tensor, 3 )
    do j = lbound( tensor, 2 ), ubound( tensor, 2 )
      write( unit, * ) tensor(:, j, i)
    end do
  end do
  close( unit )

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

  integer(i32) :: unit, i, mini, mend, nini, nend

  call assert( trim(file_format)=='text' .or. trim(file_format)=='binary', &
    'file_format must be text or binary' )
  call terminate_if_file_does_not_exist( file_name )
  select case( trim( file_format ) )
    case('text')
      call open_text_file( file_name, 'read', unit )
    case('binary')
      call open_binary_file( file_name, 'read', unit )
  end select
  read( unit, * ) i ! rank
  call terminate_if_false( i==2, 'The file ' // trim(file_name) // ' contains no matrix' )
  read( unit, * ) mini, nini, mend, nend 
  allocate( matrix(mini:mend, nini:nend) )
  do i = nini, nend 
    read( unit, * ) matrix(mini:mend, i)
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

  integer(i32) :: unit, i, j, mini, mend, nini, nend, pini, pend

  call assert( trim(file_format)=='text' .or. trim(file_format)=='binary', &
    'file_format must be text or binary' )
  call terminate_if_file_does_not_exist( file_name )
  select case( trim( file_format ) )
    case('text')
      call open_text_file( file_name, 'read', unit )
    case('binary')
      call open_binary_file( file_name, 'read', unit )
  end select
  read( unit, * ) i ! rank
  call terminate_if_false( i==3, 'The file ' // trim(file_name) // ' contains no tensor of rank 3' )
  read( unit, * ) mini, nini, pini, mend, nend, pend
  allocate( tensor(mini:mend, nini:nend, pini:pend) )
  do i = pini, pend 
    do j = nini, nend
      read( unit, * ) tensor(mini:mend, j, i)
    end do
  end do
  close( unit )

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
