!>File: matrix_contraction.f90
!> This module implements a subroutine for the second-variational style
!> contraction of complex matrices:
!>    X = factor_a * A_upper^H * B * C_upper
!>     + factor_b * A_lower^H * B * C_lower
module matrix_contraction
  use precision,  only: dp
  use constants,  only: zone, zzero
  use xlapack,    only: matrix_multiply
#include "asserts.fpp"
  implicit none

  interface contract_A_and_C_with_B
    module procedure :: contract_A_and_C_with_B_complex_dp
  end interface

contains
  !----------------------------------------------------------------------
  !> Given (A in C^{2K x M}), (B in C^{K x K}), (C in C^{2K x N}),
  !> compute:
  !>   X = factor_a * A_upper^H * B * C_upper
  !>       + factor_b * A_lower^H * B * C_lower
  !> where upper = 1:K, lower = K+1:2K
  !----------------------------------------------------------------------
  subroutine contract_A_and_C_with_B_complex_dp(A, B, C, X, factor_a, factor_b)
    !> Input arrays
    complex(dp), intent(in),  contiguous :: A(:, :), B(:, :), C(:, :)
    !> Output
    complex(dp), intent(out), contiguous :: X(:, :)
    !> Optional prefactors
    complex(dp), intent(in), optional :: factor_a, factor_b

    ! local variables
    complex(dp), allocatable :: B_times_C(:,:), A_times_B_times_C_upper(:,:), A_times_B_times_C_lower(:,:)
    complex(dp) :: factor_a_local, factor_b_local
    integer :: M, N, K

    ! shapes
    M = size(A, 2)             ! A is (2K x M)
    K = size(B, 1)             ! B is (K x K)
    N = size(C, 2)             ! C is (2K x N)

    CALL_ASSERT(size(A,1) == 2*K, "size(A,1) /= 2*K")
    CALL_ASSERT(size(B,2) == K,   "B must be KxK")
    CALL_ASSERT(size(C,1) == 2*K, "size(C,1) /= 2*K")
    CALL_ASSERT(size(X,1)==M .and. size(X,2)==N, "X must be (M x N)")

    ! set prefactors
    factor_a_local = zone
    if (present(factor_a)) factor_a_local = factor_a

    factor_b_local = zone
    if (present(factor_b)) factor_b_local = factor_b

    allocate( B_times_C(K, N), A_times_B_times_C_upper(M, N), A_times_B_times_C_lower(M, N) )

    ! upper block
    call matrix_multiply( B, C(1:K,:), B_times_C, "N", "N" )
    call matrix_multiply( A(1:K,:), B_times_C, A_times_B_times_C_upper, "C", "N" )

    ! lower block
    call matrix_multiply( B, C(K+1:2*K,:), B_times_C, "N", "N" )
    call matrix_multiply( A(K+1:2*K,:), B_times_C, A_times_B_times_C_lower, "C", "N" )

    ! combine results
    X = factor_a_local * A_times_B_times_C_upper + factor_b_local * A_times_B_times_C_lower

    deallocate(B_times_C, A_times_B_times_C_upper, A_times_B_times_C_lower)
  end subroutine contract_A_and_C_with_B_complex_dp

end module matrix_contraction
