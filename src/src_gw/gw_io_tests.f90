!> Module with unit tests for [[gw_io(module)]]
module gw_io_tests
  use gw_io, only: build_file_name, read_from_file, write_to_file, read_bounds_from_file
  use math_utils, only: all_close
  use mock_arrays, only: complex_matrix_5x7, complex_vector_7, real_vector_5
  use modmpi, only: mpiinfo
  use precision, only: i32, dp
  use unit_test_framework, only : unit_test_type

  implicit none 

  private

  character(len=*), parameter :: binary_format = 'binary'
  character(len=*), parameter :: text_format = 'text'

  public :: run_gw_io_test_driver

contains

!> Run tests for [[block_data_file(module)]]
subroutine run_gw_io_test_driver( mpiglobal, kill_on_failure )
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal
  !> Kill the program before the test driver finishes
  !> if an assertion fails
  logical, optional :: kill_on_failure
  
  type(unit_test_type) :: test_report

  call test_report%init( mpiglobal )

  ! Run and assert tests
  call test_build_file_name( test_report )
  call test_write_read( test_report, mpiglobal )
  call test_write_read_bounds( test_report, mpiglobal )

  if (present(kill_on_failure)) then
    call test_report%report( 'gw_io', kill_on_failure )
  else
    call test_report%report( 'gw_io' )
  end if

  call test_report%finalise()

end subroutine

!> Unit tests for [[build_file_name]]
subroutine test_build_file_name( test_report )
  !> Unit test report
  type(unit_test_type) :: test_report

  integer(i32), parameter :: maxlen = 60
  integer(i32), parameter :: mock_int = 4
  character(len=*), parameter :: file_extension = '.OUT'
  character(len=maxlen) :: file_name_expected
  character(len=maxlen) :: file_name
  character(len=maxlen) :: base_name

  base_name = 'test'
  call build_file_name( base_name, file_name )
  file_name_expected = trim( base_name ) // file_extension
  call test_report%assert( file_name == file_name_expected, 'String ' // trim(file_name) &
    // ' is not as expected: ' // trim(file_name_expected) )
  
  base_name = 'xyz'
  call build_file_name( base_name, mock_int, file_name )
  write( file_name_expected, * ) mock_int
  file_name_expected = trim( adjustl(base_name) ) // trim( adjustl( file_name_expected ) ) // file_extension
  call test_report%assert( file_name == file_name_expected, &
    'String ' // trim(file_name) // ' is not as expected: ' // trim(file_name_expected) )
end subroutine


!> Unit tests for [[read_from_file]] and [[write_to_file]]
subroutine test_write_read( test_report, mpiglobal )
  !> Unit test report
  type(unit_test_type) :: test_report
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal

  character(len=*), parameter :: file_extension = '.OUT'
  integer(i32), parameter :: maxlen = 50
  character(len=maxlen) :: file_name
  real(dp), parameter :: tol = 1.0e-8_dp
  real(dp), allocatable :: real_vector(:)
  complex(dp), allocatable :: vector(:), matrix(:, :), matrix_ref(:, :)
  complex(dp), allocatable :: tensor_rank3_ref(:, :, :), tensor_rank3(:, :, :)
  complex(dp), allocatable :: tensor_rank4_ref(:, :, :, :), tensor_rank4(:, :, :, :)

  ! Each MPI rank reads/writes its own file
  write( file_name, * ) mpiglobal%rank
  file_name = 'test-rank' // trim( adjustl(file_name) ) // file_extension

  ! Write/read matrix in text format
  call write_to_file( file_name, complex_matrix_5x7, lbound( complex_matrix_5x7 ), text_format )
  call read_from_file( file_name, matrix, text_format )
  call test_report%assert( all_close(matrix, complex_matrix_5x7, tol=tol), 'Matrix read from text file differs from reference' )

  ! Write/read matrix in binary format
  call write_to_file( file_name, complex_matrix_5x7, lbound( complex_matrix_5x7 ), binary_format )
  call read_from_file( file_name, matrix, binary_format )
  call test_report%assert( all_close(matrix, complex_matrix_5x7, tol=tol), 'Matrix read from binary file differs from reference' )

  ! Write/read matrix with lbounds
  allocate(matrix_ref(2:6,7:13))
  matrix_ref = complex_matrix_5x7
  call write_to_file( file_name, matrix_ref, lbound( matrix_ref ), text_format )
  call read_from_file( file_name, matrix, text_format // '  ' )
  call test_report%assert( all_close(matrix, matrix_ref, tol=tol), 'Matrix read from text file differs from reference' )
  call test_report%assert( all(lbound(matrix)==lbound(matrix_ref)), 'Matrix has wrong lbound' )

  ! Write/read matrix with lbounds (read without allocating)
  call write_to_file( file_name, matrix_ref, lbound( matrix_ref ), text_format )
  call read_from_file( file_name, matrix, lbound( matrix_ref ), text_format // ' ' )
  call test_report%assert( all_close(matrix, matrix_ref, tol=tol), 'Matrix read from text file differs from reference' )

  ! Write/read vector in binary format
  call write_to_file( file_name, complex_vector_7, binary_format )
  call read_from_file( file_name, vector, binary_format )
  call test_report%assert( all_close(vector, complex_vector_7, tol=tol), 'Vector read from binary file differs from reference' )

  ! Write/read vector in binary format (read without allocating)
  call write_to_file( file_name, complex_vector_7, 1, binary_format )
  call read_from_file( file_name, vector, 1, binary_format )
  call test_report%assert( all_close(vector, complex_vector_7, tol=tol), 'Vector read from binary file differs from reference' )

  ! Write/read real vector and complex matrix in text format
  call write_to_file( file_name, real_vector_5, matrix_ref, text_format )
  call read_from_file( file_name, real_vector, matrix, text_format )
  call test_report%assert( all_close(real_vector, real_vector_5, tol=tol), 'Vector read from text file differs from reference' )
  call test_report%assert( all_close(matrix, matrix_ref, tol=tol), 'Vector read from text file differs from reference' )

  ! Write/read tensor in binary format
  tensor_rank3_ref = reshape( complex_matrix_5x7, [2, 3, 4])
  call write_to_file( file_name, tensor_rank3_ref, lbound( tensor_rank3_ref ), binary_format )
  call read_from_file( file_name, tensor_rank3, binary_format )
  call test_report%assert( all_close(tensor_rank3, tensor_rank3_ref, tol=tol), 'Tensor read from binary file differs from reference' )
  call test_report%assert( all(lbound(tensor_rank3)==lbound(tensor_rank3_ref)), 'Tensor has wrong lbound' )

  ! Write/read tensor of rank 4 in binary format
  tensor_rank4_ref = reshape( complex_matrix_5x7, [2, 3, 2, 2])
  call write_to_file( file_name, tensor_rank4_ref, lbound( tensor_rank4_ref ), binary_format )
  call read_from_file( file_name, tensor_rank4, binary_format )
  call test_report%assert( all_close(tensor_rank4, tensor_rank4_ref, tol=tol), 'Tensor read from binary file differs from reference' )
  call test_report%assert( all(lbound(tensor_rank4)==lbound(tensor_rank4_ref)), 'Tensor has wrong lbound' )

  call delete_file( file_name )

end subroutine


!> Unit tests for [[read_bounds_from_file]] and [[write_to_file]]
subroutine test_write_read_bounds( test_report, mpiglobal )
  !> Unit test report
  type(unit_test_type) :: test_report
  !> mpi information
  type(mpiinfo), intent(in) :: mpiglobal

  character(len=*), parameter :: file_extension = '.OUT'
  integer(i32), parameter :: maxlen = 50
  character(len=maxlen) :: file_name
  complex(dp), allocatable :: vector_ref(:), matrix_ref(:, :), tensor_ref(:, :, :), tensor(:, :, :)
  real(dp), allocatable :: real_vector_ref(:)
  integer(i32) :: lbound_vector(1), ubound_vector(1)
  integer(i32) :: lbounds_matrix(2), ubounds_matrix(2)
  integer(i32) :: lbounds_tensor(3), ubounds_tensor(3)

  ! Each MPI rank reads/writes its own file
  write( file_name, * ) mpiglobal%rank
  file_name = 'test-rank' // trim( adjustl(file_name) ) // file_extension
  
  ! Write matrix in text format and read its l- and ubounds
  lbounds_matrix = [2, 7]
  ubounds_matrix = lbounds_matrix + shape(complex_matrix_5x7) - 1
  allocate(matrix_ref(lbounds_matrix(1):ubounds_matrix(1), lbounds_matrix(2):ubounds_matrix(2)))
  matrix_ref = complex_matrix_5x7
  call write_to_file( file_name, matrix_ref, lbound( matrix_ref ), text_format )
  call read_bounds_from_file( file_name, text_format, lbounds_matrix, ubounds_matrix )
  call test_report%assert( all(lbounds_matrix == lbound( matrix_ref ) ), 'Matrix lbounds do not match' )
  call test_report%assert( all(ubounds_matrix == ubound( matrix_ref ) ), 'Matrix ubounds do not match' )

  ! Write matrix in binary format and read its l- and ubounds
  call write_to_file( file_name, matrix_ref, lbound( matrix_ref ), binary_format )
  call read_bounds_from_file( file_name, binary_format, lbounds_matrix, ubounds_matrix )
  call test_report%assert( all(lbounds_matrix == lbound( matrix_ref ) ), 'Matrix lbounds do not match' )
  call test_report%assert( all(ubounds_matrix == ubound( matrix_ref ) ), 'Matrix ubounds do not match' )

  ! Write vector in binary format and read its l- and ubounds
  lbound_vector = -8
  ubound_vector = lbound_vector + shape(complex_vector_7) - 1
  allocate( vector_ref(lbound_vector(1):ubound_vector(1)) )
  vector_ref = complex_vector_7
  call write_to_file( file_name, vector_ref, lbound( vector_ref, 1 ), binary_format )
  call read_bounds_from_file( file_name, binary_format, lbound_vector, ubound_vector )
  call test_report%assert( all(lbound_vector == lbound( vector_ref ) ), 'Vector lbounds do not match' )
  call test_report%assert( all(ubound_vector == ubound( vector_ref ) ), 'Vector ubounds do not match' )

  ! Write tensor in binary format and read its l- and ubounds
  ! The size is arbitrary
  lbounds_tensor = [2, 7, 9]
  ubounds_tensor = lbounds_tensor + [5, 6, 9] - 1
  allocate(tensor_ref(lbounds_tensor(1):ubounds_tensor(1), &
                      lbounds_tensor(2):ubounds_tensor(2), &
                      lbounds_tensor(3):ubounds_tensor(3)), &
                      source=(1.0_dp, 1.0_dp) )
  call write_to_file( file_name, tensor_ref, lbound( tensor_ref ), binary_format )
  call read_bounds_from_file( file_name, binary_format, lbounds_tensor, ubounds_tensor )
  call test_report%assert( all(lbounds_tensor == lbound( tensor_ref ) ), 'Tensor lbounds do not match' )
  call test_report%assert( all(ubounds_tensor == ubound( tensor_ref ) ), 'Tensor ubounds do not match' )

  call delete_file( file_name )

end subroutine


!> Delete a file, if it exists
subroutine delete_file( file_name )
  !> File name to delete
  character(len=*), intent(in) :: file_name

  integer(i32) :: unit
  logical :: file_exists

  ! Delete test file
  inquire( file=trim(file_name), exist=file_exists )
  if( file_exists ) then
    open( newunit=unit, file=trim(file_name), status='old' )
    close( unit, status='delete' )
  end if
  
end subroutine

end module