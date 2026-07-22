!> Module with LAPACK wrappers for solving system of linear equations
!>    \[ A * x = y \]
!> where `A` is a hermitian positive definite matrix
module linear_system_positive_definite
#include "asserts.fpp"
  use lapack_f95_interfaces, only: zposv
  use math_utils, only: is_positive_definite
  use precision, only: dp, i32

  implicit none

  private

  !> Default value for uplo
  character(len=1), parameter :: uplo_default = 'U'
  character(len=1), parameter :: uplo_valid_values(4) = ['U', 'u', 'L', 'l']

  public :: positive_definite_solve

  interface positive_definite_solve
    module procedure :: positive_definite_solve_complex_dp
  end interface

contains

!> Solve the system of linear equations \(Ax=y\) for positive definite `A`
subroutine positive_definite_solve_complex_dp(A, y, uplo)
  !> Positive definite matrix `A` that determines the system of linear equation
  complex(dp), contiguous, intent(in) :: A(:, :)
  !> On entry, the right hand side matrix. On exit, the solution matrix `x`
  complex(dp), contiguous, intent(inout) :: y(:, :)
  !> Define if the part of \(A\) to be referenced:
  !> - upper triangular part: **uplo** = `'U'` or **uplo** = `'u'`
  !> - lower triangular part: **uplo** = `'L'` or **uplo** = `'l'`
  character, optional, intent(in) :: uplo

  integer(i32) :: dim_system, n_cols, info
  character :: uplo_local
  complex(dp), allocatable :: A_copy(:, :)

  if( present(uplo) ) then
    CALL_ASSERT(any(uplo == uplo_valid_values), "invalid uplo")
    uplo_local = uplo
  else
    uplo_local = uplo_default
  end if

  CALL_ASSERT( is_positive_definite(A), "A must be positive definite")
  CALL_ASSERT( size(A, 2) == size(y, 1), "A and y must have compatible sizes")

  dim_system = size(A, 1)
  n_cols = size(y, 2)
  ! Copy A because ZPOSV modifies it
  A_copy = A
  call ZPOSV( uplo_local, dim_system, n_cols, A_copy, dim_system, y, dim_system, info )
end subroutine
end module