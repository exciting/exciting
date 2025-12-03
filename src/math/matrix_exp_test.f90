!> Module for unit tests of the functions in [[matrix_exp]]
module matrix_exp_test
  use constants, only: zone, zzero, zi
  use matrix_exp
  use math_utils, only: all_close, transpose_reshape
  use modmpi, only: mpiinfo
  use precision, only: dp
  use unit_test_framework, only : unit_test_type

  implicit none
  private

  public :: matrix_exp_test_driver

  contains

  subroutine matrix_exp_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests
    call test_exp_hermitian_matrix_times_vectors( test_report )

    call test_exp_general_matrix_times_vectors( test_report )

    call test_exphouston_hermitian_matrix_times_vectors( test_report )

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('matrix_exp', kill_on_failure)
    else
      call test_report%report('matrix_exp')
    end if

    ! Finalise test object
    call test_report%finalise()

  end subroutine matrix_exp_test_driver

  !> Tests for the subroutine [[exp_hermitian_matrix_times_vectors]].
  !> 5 tests are carried out.
  subroutine test_exp_hermitian_matrix_times_vectors( test_report )
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    ! All tests are usual cases
    ! Test 01: scalar input
    call api_to_call_test( &
      & test_number=1, order_taylor=1, alpha=zone, &
      & H=reshape([zzero],[1,1]), &
      & S=reshape([zone],[1,1]), &
      & vectors=reshape([zone],[1,1]), &
      & expected=reshape([zone],[1,1]), &
      & tol=1.e-10_dp, argument_test_report=test_report )
    ! Test 02: scalar input
    call api_to_call_test( &
      & test_number=2, order_taylor=10, alpha=zone, &
      & H=reshape([zone],[1,1]), &
      & S=reshape([zone],[1,1]), &
      & vectors=reshape([zone],[1,1]), &
      & expected=reshape([2.7182818011463845_dp*zone],[1,1]), &
      & tol=1.e-10_dp, argument_test_report=test_report )
    ! Test 03: exponential of null matrix should be the identity
    call api_to_call_test( &
      & test_number=3, order_taylor=4, alpha=zone, &
      & H=reshape([zzero,zzero,zzero,zzero],[2,2]), &
      & S=reshape([zone,zzero,zzero,zone],[2,2]), &
      & vectors=reshape([zone,zzero,zzero,zone],[2,2]), &
      & expected=reshape([zone,zzero,zzero,zone],[2,2]), &
      & tol=1.e-10_dp, argument_test_report=test_report )
    ! Test 04: H, S = matrices 2x2, S = Identity
    call api_to_call_test( &
      & test_number=4, order_taylor=10, alpha=zone, &
      & H=transpose_reshape([& 
        & 0.5_dp*zone, 0.2_dp*zone,&
        & 0.2_dp*zone, 0.4_dp*zone], [2,2]), &
      & S=transpose_reshape([&
        & zone,  zzero,&
        & zzero, zone], [2,2]), &
      & vectors=transpose_reshape([&
        & zone, zzero, &
        & zzero, zone], [2,2]), &
      & expected=transpose_reshape([&
        & 1.6807292531262175_dp*zone, 0.3158889386228407_dp*zone,&
        & 0.3158889386228407_dp*zone, 1.5227847839122006_dp*zone],[2,2]), &
      & tol=1.e-10_dp, argument_test_report=test_report )
    ! Test 05: H, S = matrices 2x2, S not Identity
    call api_to_call_test( &
      & test_number=5, order_taylor=4, alpha=0.1_dp*zi, &
      & H=transpose_reshape([&
        & 0.5_dp*zone, 0.2_dp*zone, &
        & 0.2_dp*zone, 0.4_dp*zone], [2,2]), &
      & S=transpose_reshape([&
        &      zone,  -0.3_dp*zi, &
        & 0.3_dp*zi,   0.7_dp*zone], [2,2]), &
      & vectors=transpose_reshape([&
        &      zone,  zzero,&
        &     zzero,   zone], [2,2]), &
      & expected=transpose_reshape([ &
        & ( 0.987970150134193_dp, 0.056707141711965_dp ), ( -0.021040917819654_dp, 0.021694671062836_dp ), &
        & ( 0.02252502998915_dp, 0.034229059191553_dp ),  ( 1.007097835195902_dp, 0.066095152312645_dp ) ], &
        & [2,2] ), &
      & tol=1.e-10_dp, argument_test_report=test_report )
    
    contains

      !> This api is used by the subroutine
      !> [[test_exp_hermitian_matrix_times_vectors]] to avoid copy/paste.
      !> It tests if [[exp_hermitianoperator_times_wavefunctions]] provides
      !> the expected result, taking care of the error message when the result
      !> does not match the expected one.
      subroutine api_to_call_test( &
          & test_number, order_taylor, alpha, H, S, &
          & vectors, expected, tol, argument_test_report )
        !> An identifier of which test is being executed. Useful for the
        !> error message if the test does not pass.
        integer, intent(in)       :: test_number
        !> Order of the Taylor expansion used to approximate the exponential
        integer, intent(in)       :: order_taylor
        !> Complex prefactor (see [[exp_hermitianmatrix_times_vectors]])
        complex(dp), intent(in)   :: alpha
        !> Hermitian matrix
        complex(dp), intent(in)   :: H(:, :)
        !> Overlap matrix
        complex(dp), intent(in)   :: S(:, :)
        !> Array of coefficients to be tested
        complex(dp), intent(in)   :: vectors(:, :)
        !> Expected result of the test
        complex(dp), intent(in)   :: expected(:, :)
        !> Tolerance
        real(dp), intent(in)      :: tol
        !> Our test object
        type(unit_test_type), intent(inout) :: argument_test_report

        complex(dp), allocatable  :: aux(:, :)
        character(len=100)        :: error_msg

        allocate( aux(size( vectors, 1 ), size( vectors, 2 )) )
        aux = vectors
        call exp_hermitian_matrix_times_vectors( &
          & order_taylor=order_taylor, alpha=alpha, &
          & H=H, S=S, vectors=aux )
        write ( error_msg, * ) &
          & 'exp_hermitian_matrix_times_vectors, test number', &
          & test_number, ' failed.'
        call test_report%assert( all_close( a=aux , b=expected, &
          & tol=tol ), message=error_msg )
      end subroutine api_to_call_test

  end subroutine

  !> Tests for the subroutine [[exp_general_matrix_times_vectors]].
  !> 1 test is carried out.
  subroutine test_exp_general_matrix_times_vectors( test_report )
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    
    !> Vectors to test [[exp_general_matrix_times_vectors]]
    complex(dp)            :: vectors_to_test(2, 2)
    !> Expected result for the exphouston
    complex(dp)            :: expected_result(2, 2)
    !> Tolerance for comparing the expected and the obtained vectors
    real(dp), parameter    :: tol = 1.e-10_dp 

    ! Test a normal case
    vectors_to_test = transpose_reshape([&
          & -zone,        zi, &
          &  2._dp*zone,  3._dp*zone ] , [2,2] )
    expected_result = transpose_reshape([&
      & (-30.75133221112_dp, 6.97003925824_dp), (-51.65417372143_dp, 21.63493311818_dp), &
        (  4.07903080788_dp,23.05627348073_dp), ( 11.46334794687_dp, 40.19849931832_dp) ], [2,2] )
    call exp_general_matrix_times_vectors( order_taylor=4, alpha=(1._dp,1._dp), &
      H=transpose_reshape([&
        &      zone,   (1._dp,1._dp), &
        &     zzero,   0.7_dp*zone], [2,2]), &
      S=transpose_reshape([&
        &      zone,  -0.3_dp*zi, &
        & 0.3_dp*zi,   0.7_dp*zone], [2,2]), &
      vectors=vectors_to_test)
    call test_report%assert( all_close( a=vectors_to_test, &
      & b= expected_result, tol=tol ) , &
      & message='exp_general_matrix_times_vectors does not return the&
      & expected result for given input.' )

  end subroutine

  !> Tests for the subroutine [[exphouston_hermitian_matrix_times_vectors]].
  !> 1 test is carried out.
  subroutine test_exphouston_hermitian_matrix_times_vectors( test_report )
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> Vectors to test [[exphouston_hermitian_matrix_times_vectors]]
    complex(dp), allocatable :: vectors_to_test(:, :)
    !> Expected result for the exphouston
    complex(dp), allocatable :: expected_result(:, :)
    !> Tolerance for comparing the expected and the obtained vectors
    real(dp), parameter :: tol = 1.e-10_dp

    allocate( vectors_to_test(2, 2), expected_result(2, 2) )
    ! Test a case with the complete set of the eigenvectors used for expansion
    vectors_to_test = transpose_reshape([&
          & -zone,        zi, &
          &  2._dp*zone,  3._dp*zone ] , [2,2] )
    expected_result = transpose_reshape([&
      & ( -1.56610785914499_dp, 0.357804262667934_dp ), ( -0.596007992382660_dp, 1.05980026647606_dp ), &
      & ( 2.11926808755586_dp, 0.188702619714271_dp ) , ( 3.01993342215869_dp, 0.198669330793250_dp) ], [2,2] )
    call exphouston_hermitian_matrix_times_vectors( alpha=-0.1_dp*zi, &
      & H=transpose_reshape([& 
          & zone,   -zi, &
          & zi,    zone ], [2,2] ), &
      & S=transpose_reshape([& 
          &     zone,    -2._dp*zi, &
          & 2._dp*zi,   5._dp*zone ], [2,2] ), &
      & vectors=vectors_to_test, tol=tol)
    call test_report%assert( all_close( a=vectors_to_test, &
      & b= expected_result, tol=tol ) , &
      & message='exphouston_hermitian_matrix_times_vectors does not return the&
      & expected result for given input.' )
    
    ! Test a case with the partial set (1 out of 2) of the eigenvectors used for expansion
    vectors_to_test = transpose_reshape([&
          & -zone,        zi, &
          &  2._dp*zone,  3._dp*zone ] , [2,2] )
    expected_result = transpose_reshape([&
      & ( 0.83333333333862258_dp, 3.6666666666654621_dp ), ( 6.73248597155119102e-12_dp, 4.6666666666666643_dp ), &
      & ( 3.6666666666654626_dp, -0.83333333333862269_dp ), ( 4.6666666666666652_dp, -6.73248597155119183e-12_dp ) ], [2,2] )
    call exphouston_hermitian_matrix_times_vectors( alpha=-0.1_dp*zi, &
      & H=transpose_reshape([& 
          & zone,   -zi, &
          & zi,    zone ], [2,2] ), &
      & S=transpose_reshape([& 
          &     zone,    -2._dp*zi, &
          & 2._dp*zi,   4.2_dp*zone ], [2,2] ), &
      & vectors=vectors_to_test, tol=tol, n_eigs = 1)
    call test_report%assert( all_close( a=vectors_to_test, &
      & b= expected_result, tol=tol ) , &
      & message='exphouston_hermitian_matrix_times_vectors does not return the&
      & expected result for given input with n_expansion < matrix size.' )

  end subroutine

end module

