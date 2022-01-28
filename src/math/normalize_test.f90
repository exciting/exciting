!> Module for unit tests of the functions in [[matrix_exp]]
module normalize_test
  use normalize
  use math_utils, only: all_close
  use precision, only: dp
  use unit_test_framework, only : unit_test_type
  implicit none

  
  private

  public :: normalize_test_driver

  contains

  subroutine normalize_test_driver(mpiglobal, kill_on_failure)
    use modmpi, only: mpiinfo

    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 1

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_normalize_vectors( test_report )

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('normalize', kill_on_failure)
    else
      call test_report%report('normalize')
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine
  
  !> Tests for the subroutine [[normalize_vectors]].
  !> 1 test is carried out.
  subroutine test_normalize_vectors( test_report )

    use constants, only: zzero, zone, zi

    implicit none
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    complex(dp)                         :: vectors(2, 2), expected(2, 2)
    real(dp), parameter                 :: tol = 1.e-10_dp

    ! Test a normal case
    vectors = transpose( reshape([ &
      & zone,                       3._dp*zone, &
      & zi,                        -4._dp*zi], shape=[2,2] ) )
    call normalize_vectors( S=transpose(reshape([ &
      & zone,                      -2._dp*zi, &
      & 2._dp*zi,                  5._dp*zone ],shape=[2,2] )),&
      & vectors=vectors )
    expected = transpose( reshape( [&
      & 0.316227766016838_dp*zone,  0.468521285665818_dp*zone, &
      & 0.316227766016838_dp*zi,   -0.624695047554424_dp*zi ], shape=[2,2] ) )
    call test_report%assert( all_close( a=vectors, b=expected, &
      & tol=tol ), message='test_normalize_vectors failed.' )
  end subroutine test_normalize_vectors    

end module