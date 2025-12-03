!> Module for unit tests of the functions in [[matrix_exp]]
module normalize_test
  use constants, only: zzero, zone, zi
  use normalize
  use math_utils, only: all_close, transpose_reshape
  use modmpi, only: mpiinfo
  use mock_arrays, only: complex_positive_definite_matrix_5x5, complex_matrix_5x5
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
  use unit_test_framework, only : unit_test_type

  implicit none
  
  private

  public :: normalize_test_driver

  real(dp), parameter :: tol = 1.e-10_dp

  contains

  subroutine normalize_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report

    character(len=*), parameter :: module_tested = 'normalize'

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests
    call test_normalize_vectors( test_report )
    call test_norm_squared_with_positive_matrix( test_report )

    ! report results
    call test_report%report( module_tested, kill_on_failure )

    ! Finalise test object
    call test_report%finalise()
  end subroutine
  
  !> Tests for the subroutine [[normalize_vectors]].
  !> 1 test is carried out.
  subroutine test_normalize_vectors( test_report )
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    complex(dp)                         :: vectors(2, 2), expected(2, 2)

    ! Test a normal case
    vectors = transpose_reshape([ &
      & zone,                       3._dp*zone, &
      & zi,                        -4._dp*zi], [2,2] ) 
    call normalize_vectors( S=transpose_reshape([ &
      & zone,                      -2._dp*zi, &
      & 2._dp*zi,                  5._dp*zone ],[2,2] ),&
      & vectors=vectors )
    expected = transpose_reshape( [&
      & 0.316227766016838_dp*zone,  0.468521285665818_dp*zone, &
      & 0.316227766016838_dp*zi,   -0.624695047554424_dp*zi ], [2,2] ) 
    call test_report%assert( all_close( a=vectors, b=expected, &
      & tol=tol ), message='test_normalize_vectors failed.' )
  end subroutine test_normalize_vectors 
  
  !> Tests for the subroutine [[norm_squared_with_positive_matrix]].
  subroutine test_norm_squared_with_positive_matrix( test_report )
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    character(len=*), parameter :: name = 'test_norm_squared_with_positive_matrix'
    integer(i32), parameter :: m = 5, n = 3
    integer(i32) :: i, test_index
    complex(dp), parameter :: vectors(m, n) = complex_matrix_5x5(1:m, 1:n)*1e-2_dp
    complex(dp) :: tmp(m, n), identity(m, m)
    real(dp) :: norms_squared(n), expected(n)

    ! Test a generic normal case
    test_index = 1
    call norm_squared_with_positive_matrix( vectors, &
      & complex_positive_definite_matrix_5x5, norms_squared )
    tmp = matmul( complex_positive_definite_matrix_5x5, vectors )
    do i = 1, n
      expected(i) = real( dot_product( tmp(:, i), vectors(:, i) ), dp )
    end do
    call test_report%assert( all_close( norms_squared, expected, tol ), &
      name // ": test " // to_char( test_index ) // " failed." )

    ! Test a case with identity matrix
    test_index = test_index + 1
    identity = zzero
    do i = 1, m
      identity(i, i) = zone
    end do
    call norm_squared_with_positive_matrix( vectors, identity, norms_squared )
    do i = 1, n
      expected(i) = real( dot_product( vectors(:, i), vectors(:, i) ), dp )
    end do
    call test_report%assert( all_close( norms_squared, expected, tol ), &
      name // ": test " // to_char( test_index ) // " failed." )
  end subroutine 

end module

