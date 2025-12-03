!> Module for unit tests of the functions in [[projection]]
module projection_test
  use math_utils, only: all_close
  use modmpi, only: mpiinfo
  use mock_arrays, only: complex_positive_definite_matrix_5x5, complex_matrix_5x5, complex_hermitian_matrix_5x5
  use precision, only: dp, i32
  use projection, only: project_y_onto_x
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none
  
  private

  public :: projection_test_driver

  contains

  subroutine projection_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes if an assertion fails
    logical, optional, intent(in) :: kill_on_failure

    character(len=*), parameter :: module_tested = 'projection'
    type(unit_test_type) :: test_report

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests
    call test_project_vectors( test_report )

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report( module_tested, kill_on_failure )
    else
      call test_report%report( module_tested )
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine
  
  !> Tests for the subroutine [[project_y_onto_x]].
  !> 1 test is carried out.
  subroutine test_project_vectors( test_report )
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    complex(dp), allocatable :: x(:, :), y(:, :), proj(:, :), proj_ref(:, :), S(:, :), aux(:, :)
    real(dp), parameter :: tol = 1.e-9_dp
    character(len=*), parameter :: test_identifier = "test_project_vectors"
    character(len=*), parameter :: error_message = test_identifier // " failed: test number "
    integer(i32) :: m, nx, ny, test_number

    S = complex_positive_definite_matrix_5x5
    x = complex_hermitian_matrix_5x5 ! N.B. x does not necessarily need to be hermitian
    y = complex_matrix_5x5

    ! Test when all matrices have the same dimension
    test_number = 1
    nx = size(x, 2); ny = size(y, 2)
    allocate( proj(nx, ny))
    call project_y_onto_x( y, x, S, proj )
    aux = matmul(S, y)
    proj_ref = matmul( conjg(transpose(x)), aux )
    call test_report%assert( all_close( proj, proj_ref, tol ), error_message // to_char(test_number) )

    ! Test with different dims
    test_number = 2
    nx = 4; ny = 3; m = size( S, 1 )
    deallocate( proj, x, y, aux, proj_ref )
    allocate( proj(nx, ny))
    x = reshape( complex_matrix_5x5, shape=[m, nx] )
    y = reshape( complex_hermitian_matrix_5x5, shape=[m, ny] )
    call project_y_onto_x( y, x, S, proj )
    aux = matmul(S, y)
    proj_ref = matmul( conjg(transpose(x)), aux )
    call test_report%assert( all_close( proj, proj_ref, tol ), error_message // to_char(test_number) )

  end subroutine test_project_vectors    

end module