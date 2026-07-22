module projection
#include "asserts.fpp"
  use math_utils, only: is_positive_definite
  use precision, only: dp, i32
  use xlapack, only: hermitian_matrix_multiply, matrix_multiply

  implicit none

  private

  public :: project_y_onto_x

  interface project_y_onto_x
    module procedure :: project_y_onto_x_with_aux
    module procedure :: project_y_onto_x_without_aux
  end interface

contains
!> Project the vectors `y` onto `x` and store the projection coefficients as `p`.   
!> \[ p = x^\dagger S y, \]
!> where `S` is the matrix with the overlaps between basis functions.
!> \(p_{ij}\) is the projection of the `j`-th column of `y` onto the `i`-th column of `x`
subroutine project_y_onto_x_with_aux( y, x, S, proj, aux )
  !> Vectors to be projected
  complex(dp), contiguous, intent(in) :: y(:, :)
  !> Vectors onto which `y` must be projected
  complex(dp), contiguous, intent(in) :: x(:, :)
  !> Overlap matrix
  complex(dp), contiguous, intent(in) :: S(:, :)
  !> Projection coefficients of `y` onto `x`
  complex(dp), contiguous, intent(out) :: proj(:, :)
  !> Auxiliary matrix, used to compute `Sy`
  !> When `project_y_onto_x` is called inside a loop, having `aux` passed may 
  !> improve performance, as it avoids successive allocations and deallocations
  complex(dp), contiguous, intent(inout) :: aux(:, :)

  CALL_ASSERT( is_positive_definite(S), "S must be positive definite")
  associate( m => size(S, 1), n => size(y, 2) )
    CALL_ASSERT( all( shape(aux) == [m, n] ), "aux must be an m x n array")
  end associate
  call hermitian_matrix_multiply( S, y, aux )
  call matrix_multiply( x, aux, proj, trans_A='C' )
end subroutine

!> Same as `project_y_onto_x_with_aux`, but without `aux`
subroutine project_y_onto_x_without_aux( y, x, S, proj )
  !> Vectors to be projected
  complex(dp), contiguous, intent(in) :: y(:, :)
  !> Vectors onto which `y` must be projected
  complex(dp), contiguous, intent(in) :: x(:, :)
  !> Overlap matrix
  complex(dp), contiguous, intent(in) :: S(:, :)
  !> Projection coefficients of `y` onto `x`
  complex(dp), contiguous, intent(out) :: proj(:, :)

  complex(dp), allocatable :: aux(:, :)

  allocate( aux(size(S, 1), size(y, 2)) )
  call project_y_onto_x_with_aux( y, x, S, proj, aux )
end subroutine

end module