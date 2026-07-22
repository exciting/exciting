module linear_algebra_2d
  use precision, only: dp, i32

  implicit none

  private
  public :: solve_2d_cramer

  contains

  !> Solve a linear system of equations \( A \mathbf{x} = \mathbf{b} \)
  !> for a \(2 \times 2\) real matrix \( A \) using Cramer's rule.
  !> 
  !> The system is:
  !> \[
  !>   \begin{bmatrix}
  !>     a_{11} & a_{12} \\
  !>     a_{21} & a_{22}
  !>   \end{bmatrix}
  !>   \begin{bmatrix}
  !>     x_1 \\
  !>     x_2
  !>   \end{bmatrix}
  !>   =
  !>   \begin{bmatrix}
  !>     b_1 \\
  !>     b_2
  !>   \end{bmatrix}
  !> \]
  !> 
  !> The solution is obtained as:
  !> \[
  !>   \begin{aligned}
  !>   \det(A) &= a_{11} a_{22} - a_{12} a_{21}, \\
  !>   x_1 &= (b_1 a_{22} - a_{12} b_2) / \det(A), \\
  !>   x_2 &= (a_{11} b_2 - b_1 a_{21}) / \det(A).
  !>   \end{aligned}
  !> \]
  !>
!>@note
!>If the determinant is too close to zero, the routine exits with an error code info=1.
!>@endnote
  pure subroutine solve_2d_cramer(A, b, x, info)
      !> input matrix
      real(dp), intent(in) :: A(2,2)
      !> input vector
      real(dp), intent(in) :: b(2)
      !> output vector
      real(dp), intent(out) :: x(2)
      !> If info=1 det(A) is too close to 0. If info=0 everything works fine.
      integer(i32), intent(out) :: info

      ! local variables
      real(dp) :: detA, detA_inv, b_local(2)
      real(dp) :: atol

      ! relative tolerance for the determinant to be considered 0
      ! the absolute tolerance used for the comparison is calculated
      ! by multiplying rtol with the largest absolute value of A squared
      real(dp), parameter :: rtol = 1.0e-12_dp

      ! Copying, to make it save for calls like solve_2d_cramer(A, b, b, info)
      b_local(1:2) = b(1:2)

      atol = rtol * maxval(abs(A)) ** 2

      detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      if (abs(detA) <= atol) then
          info = 1
          x(1:2) = 0.0_dp
      else 
          info = 0
          detA_inv = 1.0_dp / detA
          x(1) = (b_local(1)*A(2,2) - A(1,2)*b_local(2)) * detA_inv
          x(2) = (A(1,1)*b_local(2) - b_local(1)*A(2,1)) * detA_inv
      end if

  end subroutine solve_2d_cramer


end module linear_algebra_2d