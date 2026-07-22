!> This module is designed to centralize IO operation in GW
module gw_io
#include "asserts.fpp"
  use modmpi, only: terminate_if_false
  use mod_mpi_gw, only: indexes_parallelization
  use precision, only: dp, i32, str_128
  use to_char_conversion, only: to_char

  implicit none

  private
  

  !> Number written to the header of a file to indicate that the file stores multiple arrays
  integer(i32), parameter :: indicator_file_with_multiple_arrays = -1
  !> Default extension
  character(len=*), parameter :: default_file_extension = '.OUT'
  !> Accepted file format = 'text'
  character(len=*), public, parameter :: file_format_text = 'text'
  !> Accepted file format = 'binary'
  character(len=*), public, parameter :: file_format_binary = 'binary'

  public :: build_file_name, write_to_file, read_from_file, read_bounds_from_file, &
            open_file
  
  interface write_to_file
    module procedure write_vector_to_file
    module procedure write_vector_to_file_given_lbound
    module procedure write_vector_and_matrix_to_file
    module procedure write_matrix_to_file
    module procedure write_matrix_to_file_given_lbounds
    module procedure write_tensor_of_rank_3_to_file
    module procedure write_tensor_of_rank_3_to_file_given_lbounds
    module procedure write_tensor_of_rank_4_to_file
    module procedure write_tensor_of_rank_4_to_file_given_lbounds
  end interface
  
  interface build_file_name
    module procedure build_file_name_with_integer
    module procedure build_file_name_only_adding_extension
  end interface

  interface read_from_file
    module procedure read_vector_from_file
    module procedure read_vector_from_file_without_allocating
    module procedure read_vector_and_matrix_from_file
    module procedure read_matrix_from_file
    module procedure read_matrix_from_file_without_allocating
    module procedure read_tensor_of_rank_3_from_file
    module procedure read_tensor_of_rank_4_from_file
  end interface

  interface open_file
    module procedure open_file_generic
  end interface
  
  interface write_header_to_file
    module procedure write_header_to_file_complex_array
    module procedure write_header_to_file_real_array
  end interface

contains
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


!> Same as [[write_vector_to_file_given_lbound]] with l_bound = 1
subroutine write_vector_to_file( file_name, vector, file_format )
  character(len=*), intent(in) :: file_name
  complex(dp), intent(in) :: vector(:)
  character(len=*), intent(in) :: file_format

  call write_vector_to_file_given_lbound( file_name, vector, 1, file_format )

end subroutine


!> Write a vector (array of rank 1) to a file
subroutine write_vector_to_file_given_lbound( file_name, vector, l_bound, file_format )
  !> File name where to write
  character(len=*), intent(in) :: file_name
  !> lbound of `vector`
  integer(i32), intent(in) :: l_bound
  !> Vector to be written into the file
  complex(dp), intent(in) :: vector(l_bound:)
  !> File format of output
  character(len=*), intent(in) :: file_format

  integer(i32) :: unit, i

  call open_file_generic( file_name, 'write', file_format, unit )
  call write_header_to_file( unit, file_format, vector, [l_bound] )
  select case( trim(file_format) )
    case( file_format_text )
      call print_complex_to_unit(unit, vector)
    case( file_format_binary )
      write( unit ) vector
  end select  
  close( unit )

end subroutine

!> Write a vector and a matrix to the same file
subroutine write_vector_and_matrix_to_file( file_name, vector, matrix, file_format )
  character(len=*), intent(in) :: file_name
  real(dp), intent(in) :: vector(:)
  complex(dp), intent(in) :: matrix(:, :)
  character(len=*), intent(in) :: file_format

  integer(i32), parameter :: n_arrays = 2
  integer(i32) :: unit
  
  call open_file_generic( file_name, 'write', file_format, unit )
  call write_header_file_with_multiple_arrays( unit, file_format, n_arrays )
  call write_header_to_file( unit, file_format, vector, lbounds=[1] )
  select case( trim(file_format) )
    case( file_format_text )
      write( unit, '(*(E23.17,X))' ) vector
    case( file_format_binary )
      write( unit ) vector
  end select  
  call write_contents_of_matrix_to_file_given_lbounds( unit, matrix, lbound( matrix ), file_format )
  close( unit )

end subroutine


!> Same as [[write_matrix_to_file_given_lbounds]] with lbounds = [1, 1]
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
  
  integer(i32) :: unit

  call open_file_generic( file_name, 'write', file_format, unit )
  call write_contents_of_matrix_to_file_given_lbounds( unit, matrix, lbounds, file_format ) 
  close( unit )
end subroutine


!> (private) Write a matrix (array of rank 2) to a file, given its lower bounds
!> Similar to [[write_matrix_to_file_given_lbounds]], but file must be already opened
!> and it is not closed here
subroutine write_contents_of_matrix_to_file_given_lbounds( unit, matrix, lbounds, file_format )
  !> Unit associated to the file
  integer(i32), intent(in) :: unit
  !> lbounds of `matrix`
  integer(i32), intent(in) :: lbounds(2)
  !> Matrix to be written into the file
  complex(dp), intent(in) :: matrix(lbounds(1):, lbounds(2):)
  !> File format of output
  character(len=*), intent(in) :: file_format
  
  integer(i32) :: i

  call write_header_to_file( unit, file_format, matrix, lbounds )
  select case( trim(file_format) )
    case( file_format_text )
      do i = lbounds(2), ubound( matrix, 2 )
        call print_complex_to_unit(unit, matrix(:, i))
      end do
    case( file_format_binary )
      do i = lbounds(2), ubound( matrix, 2 )
        write( unit ) matrix(:, i)
      end do
  end select  
end subroutine
  

!> Same as [[write_tensor_of_rank_3_to_file_given_lbounds]] with lbounds = [1, 1, 1]
subroutine write_tensor_of_rank_3_to_file( file_name, tensor, file_format )
  !> File name where to write
  character(len=*), intent(in) :: file_name
  !> Tensor to be written into the file
  complex(dp), intent(in) :: tensor(:, :, :)
  !> Format of the output file
  character(len=*), intent(in) :: file_format

  call write_tensor_of_rank_3_to_file_given_lbounds( file_name, tensor, [1, 1, 1], file_format )

end subroutine


!> Write a tensor of rank 3 to a file
subroutine write_tensor_of_rank_3_to_file_given_lbounds( file_name, tensor, lbounds, file_format )
  character(len=*), intent(in) :: file_name
  integer(i32), intent(in) :: lbounds(3)
  complex(dp), intent(in) :: tensor(lbounds(1):, lbounds(2):, lbounds(3):)
  character(len=*), intent(in) :: file_format

  integer(i32) :: unit, i, j

  call open_file_generic( file_name, 'write', file_format, unit )
  call write_header_to_file( unit, file_format, tensor, lbounds )
  select case( trim(file_format) )
    case( file_format_text )
      do i = lbounds(3), ubound( tensor, 3 )
        do j = lbounds(2), ubound( tensor, 2 )
          call print_complex_to_unit(unit, tensor(:, j, i))
        end do
      end do
    case( file_format_binary )
      do i = lbounds(3), ubound( tensor, 3 )
        do j = lbounds(2), ubound( tensor, 2 )
          write( unit ) tensor(:, j, i)
        end do
      end do
  end select
  close( unit )

end subroutine


!> Same as [[write_tensor_of_rank_4_to_file_given_lbounds]] with lbounds = [1, 1, 1, 1]
subroutine write_tensor_of_rank_4_to_file( file_name, tensor, file_format )
  !> File name where to write
  character(len=*), intent(in) :: file_name
  !> Tensor to be written into the file
  complex(dp), intent(in) :: tensor(:, :, :, :)
  !> Format of the output file
  character(len=*), intent(in) :: file_format

  call write_tensor_of_rank_4_to_file_given_lbounds( file_name, tensor, [1, 1, 1, 1], file_format )

end subroutine


!> Write a tensor of rank 4 to a file
subroutine write_tensor_of_rank_4_to_file_given_lbounds( file_name, tensor, lbounds, file_format )
  character(len=*), intent(in) :: file_name
  integer(i32), intent(in) :: lbounds(4)
  complex(dp), intent(in) :: tensor(lbounds(1):, lbounds(2):, lbounds(3):, lbounds(4):)
  character(len=*), intent(in) :: file_format

  integer(i32) :: unit, i, j, k

  call open_file_generic( file_name, 'write', file_format, unit )
  call write_header_to_file( unit, file_format, tensor, lbounds )
  select case( trim(file_format) )
    case( file_format_text )
      do i = lbounds(4), ubound( tensor, 4 )
        do j = lbounds(3), ubound( tensor, 3 )
          do k = lbounds(2), ubound( tensor, 2 )
            call print_complex_to_unit(unit, tensor(:, k, j, i))
          end do
        end do
      end do
    case( file_format_binary )
      do i = lbounds(4), ubound( tensor, 4 )
        do j = lbounds(3), ubound( tensor, 3 )
          do k = lbounds(2), ubound( tensor, 2 )
            write( unit ) tensor(:, k, j, i)
          end do
        end do
      end do
  end select
  close( unit )

end subroutine

!> (private) Write a header to an output file
!> to indicate that this file stores multiple arrays
subroutine write_header_file_with_multiple_arrays( unit, file_format, n_arrays )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Format of the file
  character(len=*), intent(in) :: file_format
  !> Number of arrays to be stored in the outpur file
  integer(i32), intent(in) :: n_arrays

  select case( trim(file_format) )
    case( file_format_text )
      write( unit, * ) indicator_file_with_multiple_arrays
      write( unit, * ) n_arrays
    case( file_format_binary )
      write( unit ) indicator_file_with_multiple_arrays
      write( unit ) n_arrays
  end select
end subroutine


!> (private) Write a header to an output file
!> The header contains the rank of the array, its lbounds and ubounds
subroutine write_header_to_file_complex_array( unit, file_format, array, lbounds )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Format of the file
  character(len=*), intent(in) :: file_format
  !> Array to be written into the file
  complex(dp), intent(in) :: array(..)
  !> Lbounds of `array`
  integer(i32), intent(in) :: lbounds(:)

  ! Note for developers:
  ! Some compilers pass the full decriptor 
  ! for assumed rank arrays. This means
  ! that the bounds are preserved. We 
  ! do a check for these cases to
  ! properly save the bounds.
  !
  ! This means that ubound( array ) + lbounds - 1
  ! is not safe

  CALL_ASSERT( size(lbounds) == rank(array), 'Incompatible array rank and lbounds' )
  select case( trim(file_format) )
    case( file_format_text )
      write( unit, * ) rank( array )
      write( unit, * ) lbounds, shape( array ) + lbounds - 1
    case( file_format_binary )
      write( unit ) rank( array )
      write( unit ) lbounds, shape( array ) + lbounds - 1
  end select

end subroutine


!> (private) Write a header to an output file
!> The header contains the rank of the array, its lbounds and ubounds
subroutine write_header_to_file_real_array( unit, file_format, array, lbounds )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Format of the file
  character(len=*), intent(in) :: file_format
  !> Array to be written into the file
  real(dp), intent(in) :: array(..)
  !> Lbounds of `array`
  integer(i32), intent(in) :: lbounds(:)

  CALL_ASSERT( size(lbounds) == rank(array), 'Incompatible array rank and lbounds' )
  select case( trim(file_format) )
    case( file_format_text )
      write( unit, * ) rank( array )
      write( unit, * ) lbounds, shape( array ) + lbounds - 1
    case( file_format_binary )
      write( unit ) rank( array )
      write( unit ) lbounds, shape( array ) + lbounds - 1
  end select

end subroutine


!> (private) Read the header of an input file
!> The header contains the rank of the array, its lbounds and ubounds
subroutine read_header_of_file( unit, file_format, rank_of_array, lbounds, ubounds )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Format of the file
  character(len=*), intent(in) :: file_format
  !> Rank of the array stored in the file
  integer(i32), intent(out) :: rank_of_array
  !> Lbounds of the array
  integer(i32), intent(out) :: lbounds(:)
  !> Ubounds of the array
  integer(i32), intent(out) :: ubounds(:)
  
  CALL_ASSERT( size(lbounds) == size(ubounds), 'lbounds and ubounds must have same size')
  select case( trim(file_format) )
    case( file_format_text )
      read( unit, * ) rank_of_array
      call terminate_if_false( rank_of_array==size(lbounds), &
        'Incompatible rank of array. Stored in file is: ' // to_char(rank_of_array) // '; required for this calculation is: ' // to_char( size(lbounds) )  )
      read( unit, * ) lbounds, ubounds
    case( file_format_binary )
      read( unit ) rank_of_array
      call terminate_if_false( rank_of_array==size(lbounds), &
        'Incompatible rank of array. Stored in file is: ' // to_char(rank_of_array) // '; required for this calculation is: ' // to_char( size(lbounds) )  )
      read( unit ) lbounds, ubounds
  end select

end subroutine


!> (private) Read a header from a file that this file stores multiple arrays
subroutine read_header_file_with_multiple_arrays( unit, file_format, n_arrays )
  !> Unit associated to the file where to write
  integer(i32), intent(in) :: unit
  !> Format of the file
  character(len=*), intent(in) :: file_format
  !> Number of arrays to be stored in the outpur file
  integer(i32), intent(out) :: n_arrays

  integer(i32) :: indicator

  select case( trim(file_format) )
    case( file_format_text )
      read( unit, * ) indicator
      read( unit, * ) n_arrays
    case( file_format_binary )
      read( unit ) indicator
      read( unit ) n_arrays
  end select

  call terminate_if_false( indicator==indicator_file_with_multiple_arrays, 'File does not contain multiple arrays')
end subroutine

!> Check if a file exists and terminates the execution if it does not
subroutine terminate_if_file_does_not_exist( file_name )
  !> Name of the file to check
  character(len=*), intent(in) :: file_name

  logical :: ok
  inquire( file=trim(file_name), exist=ok )
  call terminate_if_false( ok, 'File '//trim(file_name)//' not found')
end subroutine
  
!> Read a complex vector (array of rank 1) from a file
subroutine read_vector_from_file( file_name, vector, file_format )
  !> Name of the file where the vector is stored
  character(len=*), intent(in) :: file_name
  !> Vector where to save the data read from file
  complex(dp), intent(out), allocatable :: vector(:)
  !> Format of the file
  character(len=*), intent(in) :: file_format

  integer(i32) :: unit, read_rank, lbound_(1), ubound_(1)
  integer(i32), parameter :: expected_rank = 1 !rank of a vector

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, file_format, read_rank, lbound_, ubound_ )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no vector' )
  allocate( vector(lbound_(1):ubound_(1)) )
  call read_vector_from_file_header_skipped( unit, vector, file_format )

end subroutine


!> Read a vector (array of rank 1) from a file and store the data in an already allocated vector
subroutine read_vector_from_file_without_allocating( file_name, vector, l_bound, file_format )
  !> Name of the file where the matrix is stored
  character(len=*), intent(in)          :: file_name
  !> Lower bounds of the matrix
  integer(i32), intent(in)              :: l_bound
  !> Vector where to save the data read from file
  complex(dp), intent(out)              :: vector(l_bound:)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: unit, l_bound_(1), u_bound_(1), read_rank
  integer(i32), parameter :: expected_rank = 1 !rank of a vector

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, file_format, read_rank, l_bound_, u_bound_ )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no vector' )
  call terminate_if_false( l_bound_(1)==l_bound .and. u_bound_(1)==ubound(vector,1) , &
    'The file ' // trim(file_name) // ' has invalid bounds' )
  call read_vector_from_file_header_skipped( unit, vector, file_format )

end subroutine


!> (private) Read a vector (array of rank 1) from a file
!> It is assumed that the file has been opened and the header has already been read
subroutine read_vector_from_file_header_skipped( unit, vector, file_format )
  !> Unit corresponding to the file where the vector is stored
  integer(i32), intent(in) :: unit
  !> Vector where to save the data read from file
  complex(dp), intent(out) :: vector(:)
  !> Format of the file
  character(len=*), intent(in) :: file_format

  select case( trim(file_format) )
    case( file_format_text )
      read( unit, * ) vector
    case( file_format_binary )
      read( unit ) vector
  end select 
  close( unit )

end subroutine


!> Read a vector and a matrix from the same file, which must have been written according to [[write_vector_and_matrix_to_file]]
subroutine read_vector_and_matrix_from_file( file_name, vector, matrix, file_format )
  !> Name of the file where the matrix is stored
  character(len=*), intent(in) :: file_name
  !> Vector where to save the data read from file
  real(dp), intent(out), allocatable :: vector(:)
  !> Matrix to store the data read from the file
  complex(dp), intent(out), allocatable :: matrix(:, :)
  !> Format of the file
  character(len=*), intent(in) :: file_format

  integer(i32), parameter :: n_arrays_expected = 2
  integer(i32), parameter :: rank_expected_vector = 1
  integer(i32), parameter :: rank_expected_matrix = 2
  integer(i32) :: unit, n_arrays, rank_read
  integer(i32) :: lbounds_vector(1), ubounds_vector(1), lbounds_matrix(2), ubounds_matrix(2)
  
  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_file_with_multiple_arrays( unit, file_format, n_arrays )
  call terminate_if_false( n_arrays == n_arrays_expected, 'File ' // trim( file_name ) // &
    'does not contain ' // to_char(n_arrays_expected) // ' arrays' )
  call read_header_of_file( unit, file_format, rank_read, lbounds_vector, ubounds_vector )
  call terminate_if_false( rank_read == rank_expected_vector, 'The file ' // trim(file_name) // ' contains no vector' )
  allocate( vector(lbounds_vector(1):ubounds_vector(1)) )
  select case( trim(file_format) )
    case( file_format_text )
      read( unit, * ) vector
    case( file_format_binary )
      read( unit ) vector
  end select  
  call read_matrix_from_already_opened_file( unit, matrix, file_format )

end subroutine


!> Read a matrix (array of rank 2) from a file
subroutine read_matrix_from_file( file_name, matrix, file_format )
  !> Name of the file to read
  character(len=*), intent(in)          :: file_name
  !> Matrix to store the data read from the file
  complex(dp), intent(out), allocatable :: matrix(:, :)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: unit

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_matrix_from_already_opened_file( unit, matrix, file_format )

end subroutine

!> (private) Read a matrix (array of rank 2) from an already opened file
subroutine read_matrix_from_already_opened_file( unit, matrix, file_format )
  !> Unit corresponding to the file where the matrix is stored
  integer(i32), intent(in) :: unit
  !> Matrix to store the data read from the file
  complex(dp), intent(out), allocatable :: matrix(:, :)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: read_rank, lbounds(2), ubounds(2)
  integer(i32), parameter :: expected_rank = 2 !rank of a matrix

  call read_header_of_file( unit, file_format, read_rank, lbounds, ubounds )
  call terminate_if_false( read_rank==expected_rank, 'The file with unit' // to_char( unit ) // ' contains no matrix' )
  allocate( matrix(lbounds(1):ubounds(1), lbounds(2):ubounds(2)) )
  call read_matrix_from_file_header_skipped( unit, matrix, file_format )

end subroutine


!> Read a matrix (array of rank 2) from a file and store the data in an already allocated matrix
subroutine read_matrix_from_file_without_allocating( file_name, matrix, lbounds, file_format )
  !> Name of the file where the matrix is stored
  character(len=*), intent(in)          :: file_name
  !> Lower bounds of the matrix
  integer(i32), intent(in)              :: lbounds (2)
  !> Vector where to save the data read from file
  complex(dp), intent(out)              :: matrix(lbounds(1):, lbounds(2):)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: unit, read_rank, lbounds_(2), ubounds_(2)
  integer(i32), parameter :: expected_rank = 2 !rank of a matrix

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, file_format, read_rank, lbounds_, ubounds_ )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no matrix' )
  call terminate_if_false( all(lbounds_==lbounds) .and. all( ubounds_ == ubound(matrix) ), &
    'The file ' // trim(file_name) // ' has invalid bounds' )
  call read_matrix_from_file_header_skipped( unit, matrix, file_format )

end subroutine


!> (private) Read a matrix (array of rank 2) from a file
!> It is assumed that the file has been opened and the header has already been read
subroutine read_matrix_from_file_header_skipped( unit, matrix, file_format )
  !> Unit corresponding to the file where the matrix is stored
  integer(i32), intent(in) :: unit
  !> Matrix where to save the data read from file
  complex(dp), intent(out) :: matrix(:, :)
  !> Format of the file
  character(len=*), intent(in) :: file_format

  integer(i32) :: i

  select case( trim(file_format) )
    case( file_format_text )
      do i = 1, size( matrix, 2 )
        read( unit, * ) matrix(:, i)
      end do
    case( file_format_binary )
      do i = 1, size( matrix, 2 )
        read( unit ) matrix(:, i)
      end do
  end select 
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
  call read_header_of_file( unit, file_format, read_rank, lbounds, ubounds )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no tensor of rank 3' )
  allocate( tensor(lbounds(1):ubounds(1), lbounds(2):ubounds(2), lbounds(3):ubounds(3)) )
  select case( trim(file_format) )
    case( file_format_text )
      do i = lbounds(3), ubounds(3) 
        do j = lbounds(2), ubounds(2)
          read( unit, * ) tensor(:, j, i)
        end do
      end do
    case( file_format_binary )
      do i = lbounds(3), ubounds(3) 
        do j = lbounds(2), ubounds(2)
          read( unit ) tensor(:, j, i)
        end do
      end do
  end select 
  close( unit )

end subroutine


!> Read a tensor of rank 4 stored in a file
subroutine read_tensor_of_rank_4_from_file( file_name, tensor, file_format )
  !> Name of the file to read
  character(len=*), intent(in)          :: file_name
  !> Tensor to store the data read from the file
  complex(dp), intent(out), allocatable :: tensor(:, :, :, :)
  !> Format of the file
  character(len=*), intent(in)          :: file_format

  integer(i32) :: unit, i, j, k, read_rank, lbounds(4), ubounds(4)
  integer(i32), parameter :: expected_rank = 4 !rank of tensor

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, file_format, read_rank, lbounds, ubounds )
  call terminate_if_false( read_rank==expected_rank, 'The file ' // trim(file_name) // ' contains no tensor of rank 4' )
  allocate( tensor(lbounds(1):ubounds(1), lbounds(2):ubounds(2), lbounds(3):ubounds(3), lbounds(4):ubounds(4)) )
  select case( trim(file_format) )
    case( file_format_text )
      do i = lbounds(4), ubound( tensor, 4 )
        do j = lbounds(3), ubound( tensor, 3 )
          do k = lbounds(2), ubound( tensor, 2 )
            read( unit, * ) tensor(:, k, j, i)
          end do
        end do
      end do
    case( file_format_binary )
      do i = lbounds(4), ubound( tensor, 4 )
        do j = lbounds(3), ubound( tensor, 3 )
          do k = lbounds(2), ubound( tensor, 2 )
            read( unit ) tensor(:, k, j, i)
          end do
        end do
      end do
  end select 
  close( unit )

end subroutine


!> Read the lower and upper bounds of an array from a file
subroutine read_bounds_from_file( file_name, file_format, lbounds, ubounds )
  !> Name of the file where the matrix is stored
  character(len=*), intent(in)          :: file_name
  !> Format of the file
  character(len=*), intent(in)          :: file_format
  !> Lower bounds of the matrix
  integer(i32), intent(out)              :: lbounds(:)
  !> Upper bounds of the matrix
  integer(i32), intent(out)              :: ubounds(:)

  integer(i32) :: unit, rank_of_array

  call open_file_generic( file_name, 'read', file_format, unit )
  call read_header_of_file( unit, file_format, rank_of_array, lbounds, ubounds )
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

  CALL_ASSERT( trim(file_format)==file_format_text .or. trim(file_format)==file_format_binary,  'file_format must be '// file_format_text // ' or ' // file_format_binary )
  
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

  open( newunit=unit, file=trim(file_name), action=action, form='formatted' )

end subroutine

  
!> Open a binary file
subroutine open_binary_file( file_name, action, unit )
  !> Name of the file
  character(len=*), intent(in)  :: file_name
  !> Action: should be e. g. "read", "write"
  character(len=*), intent(in)  :: action
  !> Unit of the file (that will be opened)
  integer(i32), intent(out)     :: unit

  open( newunit=unit, file=trim(file_name), action=action, form='unformatted', access='stream' )
  
end subroutine

subroutine print_complex_to_unit(unit, a)
    integer(i32), intent(in) :: unit
    complex(dp), intent(in) :: a(:)

    integer(i32) :: i

    do i = 1, size(a)
        write(unit,'("(",E23.17,",",E23.17,")",X)', advance="no") real(a(i)), aimag(a(i))
    end do
    write(unit,'()', advance="yes")

end subroutine print_complex_to_unit
  
end module
