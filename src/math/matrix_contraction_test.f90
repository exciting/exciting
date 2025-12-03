
module matrix_contraction_test
  use precision,            only: dp
  use modmpi,               only: mpiinfo
  use unit_test_framework,  only: unit_test_type
  use matrix_contraction,      only: contract_A_and_C_with_B
  use xlapack,              only: matrix_multiply
  use math_utils,           only: fill_random,all_close
  implicit none

  private
  public :: matrix_contraction_test_driver

contains

!> Driver for [[matrix_contractionion]] unit tests.
  subroutine matrix_contraction_test_driver(mpiglobal, kill_on_failure)
    type(mpiinfo), intent(in) :: mpiglobal
    logical, intent(in), optional :: kill_on_failure

    type(unit_test_type) :: test_report

    call test_report%init( mpiglobal)

    call test_small_case(test_report)
    call test_random_case(test_report)
    call test_random_case_with_factors(test_report)

    ! Report
    if (present(kill_on_failure)) then
       call test_report%report("matrix_contraction", kill_on_failure)
    else
       call test_report%report("matrix_contraction")
    end if

    call test_report%finalise()
  end subroutine matrix_contraction_test_driver

  subroutine test_small_case(test)
    class(unit_test_type), intent(inout) :: test

    integer, parameter :: k=1, m=2, n=2
    complex(dp) :: A(2*k,m), B(k,k), C(2*k,n)
    complex(dp) :: X(m,n), X_ref(m,n)

    ! Prepare data
    A(:,:) = (1.0_dp, 0.0_dp)
    B(1,1) = (2.0_dp, 0.0_dp)
    C(:,:) = (0.0_dp, 0.0_dp)
    C(1,:) = (3.0_dp, 0.0_dp)

    X_ref(:,:) = (6.0_dp, 0.0_dp)

    call contract_A_and_C_with_B(A, B, C, X)

    call test%assert( all_close(X, X_ref), &
         "test_small_case: X not matching expected 6.0" )
  end subroutine test_small_case

 !> A random test without explicit factors (k=4,m=5,n=6).
  subroutine test_random_case(test)
    class(unit_test_type), intent(inout) :: test

    integer, parameter :: k=4, m=5, n=6
    complex(dp) :: A(2*k,m), B(k,k), C(2*k,n)
    complex(dp) :: X(m,n), X_manual(m,n)
    complex(dp) :: B_times_C_upper(k, n), X_upper(m, n)
    complex(dp) :: B_times_C_lower(k, n), X_lower(m, n)

    call fill_random(A)
    call fill_random(B)
    call fill_random(C)

 
    call contract_A_and_C_with_B(A, B, C, X)

    B_times_C_upper = matmul(B, C(1:k, :))
    X_upper         = matmul(conjg(transpose(A(1:k, :))), B_times_C_upper)

    B_times_C_lower = matmul(B, C(k+1:2*k, :))
    X_lower         = matmul(conjg(transpose(A(k+1:2*k, :))), B_times_C_lower)

    X_manual = X_upper + X_lower

    call test%assert(all_close(X, X_manual, 1.0d-12), &
         "test_random_case: mismatch between code and manual approach")
  end subroutine test_random_case

  !> A random test WITH factor_a, factor_b
   subroutine test_random_case_with_factors(test)
    class(unit_test_type), intent(inout) :: test

    integer, parameter :: k=3, m=4, n=5
    complex(dp), parameter :: factor_a = cmplx(1.2_dp, -0.3_dp)
    complex(dp), parameter :: factor_b = cmplx(0.7_dp,  0.9_dp)

    complex(dp) :: A(2*k,m), B(k,k), C(2*k,n)
    complex(dp) :: X(m,n), X_manual(m,n)
    complex(dp) :: B_times_C_upper(k,n), B_times_C_lower(k,n)
    complex(dp) :: X_upper(m,n),  X_lower(m,n)

    call fill_random(A)
    call fill_random(B)
    call fill_random(C)

    call contract_A_and_C_with_B(A, B, C, X, factor_a, factor_b)
    !    upper block
    B_times_C_upper = matmul(B, C(1:k, :))
    X_upper         = matmul(conjg(transpose(A(1:k, :))), B_times_C_upper)
    X_upper         = factor_a * X_upper    
    !    lower block
    B_times_C_lower = matmul(B, C(k+1:2*k, :))
    X_lower         = matmul(conjg(transpose(A(k+1:2*k, :))), B_times_C_lower)
    X_lower         = factor_b * X_lower    
    
    X_manual = X_upper + X_lower

    call test%assert(all_close(X, X_manual, 1.0d-12), &
         "test_random_case_with_factors: mismatch with factors")
  end subroutine test_random_case_with_factors

end module matrix_contraction_test
