!> Module for unit tests of the functions in [[advanced_matrix_operations]]
module integration_test
  use integration
  use math_utils, only: all_close
  use precision, only: dp
  use unit_test_framework, only : unit_test_type

  implicit none
  private

  public :: integration_test_driver

  contains

  subroutine integration_test_driver(mpiglobal, kill_on_failure)
    use modmpi, only: mpiinfo

    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 2

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_ODESolver_RungeKutta4thOrder( test_report )

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('integration', kill_on_failure)
    else
      call test_report%report('integration')
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine


  !> Tests for the subroutine [[ODESolver_RungeKutta4thOrder]].
  !> 2 tests are carried out.
  subroutine test_ODESolver_RungeKutta4thOrder( test_report )
    use constants, only: zi, zone, zzero

    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    ! Wavefunction to be evolved by [[ODESolver_RungeKutta4thOrder]]
    complex(dp)                         :: psi(3, 2)
    ! Expected result after calling [[ODESolver_RungeKutta4thOrder]]
    complex(dp)                         :: expected_result(3 ,2)
    ! Auxiliary matrices
    complex(dp)                         :: H(3,3), H_future(3,3)
    ! Tolerance for comparing the expected and the obtained vectors
    real(dp), parameter                 :: tol = 1.e-10_dp
    

    ! Test a normal case
    psi = transpose(reshape([&
          & 9._dp*zone,    0.7_dp*zone, &
          &    3_dp*zi,    -2._dp*zone, &
          &      zzero,           zone], [2, 3] ))
    expected_result = transpose(reshape([&
      & (9.14980501921296_dp, -0.372522280864198_dp),     (0.710213643595679_dp, -0.290095949485597_dp), &
      & (-8.105717901234555e-2_dp, 3.46144749591049_dp ), (-1.87477479948560_dp, 0.491909162309671_dp), &
      & zzero,                                            (0.870316026041667_dp, -0.492234854166667_dp) ], &
      & [2, 3] ))
    call ODESolver_RungeKutta4thOrder( time_step=0.1_dp, alpha=zi, &
      & H=transpose(reshape([&
          & zone,           zi,        zzero, &
          &  -zi,   2._dp*zone,        zzero, &
          & zzero,       zzero,   5._dp*zone ], [3, 3] )), &
      & H_past=transpose(reshape([& 
          &  0.9_dp*zone,      1.1_dp*zi,            zzero,  &
          &   -1.1_dp*zi,    2.1_dp*zone,            zzero,  &
          &        zzero,          zzero,      4.7_dp*zone ], [3, 3] )), &
      & S=transpose(reshape([& 
          &     4._dp*zone,   (2._dp,1_dp),    zzero, &
          & (2._dp,-1._dp),      2_dp*zone,    zzero, &
          &          zzero,          zzero,     zone], [3, 3] )), &
      & x=psi )
    call test_report%assert( all_close( a=psi, b=expected_result, tol=tol ), &
      & message='test_ODESolver_RungeKutta4thOrder - First test failed.')

    ! Test when making we obtain by hand H(t-dt) from H(t+dt) and H(t)
    psi = transpose(reshape([&
        & 9._dp*zone,   0.7_dp*zone, &
        &    3_dp*zi,   -2._dp*zone, &
        &      zzero,          zone], [2, 3] ))
    expected_result=transpose(reshape([&
      & (9.19811194598765_dp, -0.278502205555555_dp),  (0.707046930002572_dp, -0.296011989670782_dp), &
      & (-0.210419102469136_dp, 3.39670194104938_dp),  (-1.87012088066872_dp, 0.506675206342592_dp) , &
      &  zzero,                                        (0.884690744791667_dp, -0.465985979166667_dp)], [2, 3] ))
    H=transpose(reshape([& 
      &  0.9_dp*zone,      1.1_dp*zi,            zzero,  &
      &   -1.1_dp*zi,    2.1_dp*zone,            zzero,  &
      &        zzero,          zzero,      4.7_dp*zone ], [3, 3] ))
    H_future = transpose(reshape([& 
      &  zone,           zi,        zzero, &
      &   -zi,   2._dp*zone,        zzero, &
      & zzero,        zzero,   5._dp*zone ], [3, 3] ))
    call ODESolver_RungeKutta4thOrder( time_step=0.1_dp, alpha=zi, &
      & H_past=2*H-H_future, &
      & H=H, &
      & S=transpose(reshape([& 
          &    4._dp*zone,       (2._dp,1_dp),     zzero, &
          & (2._dp,-1._dp),         2_dp*zone,     zzero, &
          &         zzero,              zzero,      zone], [3, 3] )), &
      & x=psi )
    call test_report%assert( all_close( a=psi, b=expected_result, tol=tol ), &
      message='test_ODESolver_RungeKutta4thOrder - Second test failed.')
  end subroutine

end module integration_test