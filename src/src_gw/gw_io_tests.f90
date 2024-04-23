!> Module with unit tests for [[gw_io(module)]]
module gw_io_tests
  use gw_io, only: build_file_name, read_from_file, write_to_file
  use math_utils, only: all_close
  use mock_arrays, only: complex_matrix_5x7, complex_vector_7
  use modmpi, only: mpiinfo
  use precision, only: i32, dp
  use unit_test_framework, only : unit_test_type

  implicit none 

  private

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
  integer, parameter :: n_assertions = 2+7

  call test_report%init( n_assertions, mpiglobal )

  ! Run and assert tests
  call test_build_file_name( test_report )
  call test_write_read( test_report, mpiglobal )

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
  integer(i32) :: unit
  logical :: file_exists
  complex(dp), allocatable :: vector(:), matrix(:, :), matrix_ref(:, :)
  complex(dp), allocatable :: tensor_ref(:, :, :), tensor(:, :, :)

  ! Each MPI rank reads/writes its own file
  write( file_name, * ) mpiglobal%rank
  file_name = 'test-rank' // trim( adjustl(file_name) ) // file_extension

  ! Write/read matrix in text format
  call write_to_file( file_name, complex_matrix_5x7, lbound( complex_matrix_5x7 ), 'text' )
  call read_from_file( file_name, matrix, 'text' )
  call test_report%assert( all_close(matrix, complex_matrix_5x7, tol=tol), 'Matrix read from text file differs from reference' )

  ! Write/read matrix in binary format
  call write_to_file( file_name, complex_matrix_5x7, lbound( complex_matrix_5x7 ), 'binary' )
  call read_from_file( file_name, matrix, 'binary' )
  call test_report%assert( all( matrix == complex_matrix_5x7 ) , 'Matrix read from binary file differs from reference' )

  ! Write/read matrix with lbounds
  allocate(matrix_ref(2:6,7:13))
  matrix_ref = complex_matrix_5x7
  call write_to_file( file_name, matrix_ref, lbound( matrix_ref ), 'text' )
  call read_from_file( file_name, matrix, 'text  ' )
  call test_report%assert( all_close(matrix, matrix_ref, tol=tol), 'Matrix read from text file differs from reference' )
  call test_report%assert( all(lbound(matrix)==lbound(matrix_ref)), 'Matrix has wrong lbound' )

  ! Write/read vector in binary format
  call write_to_file( file_name, complex_vector_7, 'binary' )
  call read_from_file( file_name, vector, 'binary' )
  call test_report%assert( all( vector == complex_vector_7 ), 'Vector read from binary file differs from reference' )

  ! Write/read tensor in binary format
  tensor_ref = reshape( complex_matrix_5x7, [2, 3, 4])
  call write_to_file( file_name, tensor_ref, lbound( tensor_ref ), 'binary' )
  call read_from_file( file_name, tensor, 'binary' )
  call test_report%assert( all( tensor == tensor_ref ), 'Tensor read from binary file differs from reference' )
  call test_report%assert( all(lbound(tensor)==lbound(tensor_ref)), 'Tensor has wrong lbound' )

  ! Delete test file
  inquire( file=file_name, exist=file_exists )
  if( file_exists ) then
    open( newunit=unit, file=file_name, status='old' )
    close( unit, status='delete' )
  end if
end subroutine

end module