module linear_algebra_3d
  use asserts, only: assert
  use precision, only: dp
  use constants, only: pi
  use math_utils, only: all_zero

  implicit none

  private
  public :: determinant_3d, &
            cross_product, &
            inverse_3d, &
            triple_product

  !> Calculate the determinant of a matrix \( \mathbf{A} \).
  interface determinant_3d
    module procedure :: determinant_3d_real_dp, determinant_3d_integer
  end interface

  !> Calculate triple product for three vectors \( \mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3 \in \mathbb{R}^3 \):
  !> \[
  !>   s = \mathbf{a}_{1} \cdot (\mathbf{a}_{2} \times \mathbf{a}_{3}).
  !> \]
  !> The routine takes the input vectors independently or columnwise as \(3 \cross 3\) matrix.
  interface triple_product
    module procedure :: determinant_3d_real_dp, triple_product_vector_inputs
  end interface

  contains 

  ! Calculate the determinant of a matrix \( \mathbf{A} \).
  pure real(dp) function determinant_3d_real_dp(A)
    !> Input matrix
    real(dp), intent(in) :: A(3, 3)

    determinant_3d_real_dp = A(1, 1) * A(2, 2) * A(3, 3)  -  A(3, 1) * A(2, 2) * A(1, 3) &
                           + A(1, 2) * A(2, 3) * A(3, 1)  -  A(3, 2) * A(2, 3) * A(1, 1) &
                           + A(1, 3) * A(2, 1) * A(3, 2)  -  A(3, 3) * A(2, 1) * A(1, 2)
  end function determinant_3d_real_dp

  ! Calculate the determinant of a matrix \( \mathbf{A} \).
  pure integer function determinant_3d_integer(A)
    !> Input matrix
    integer, intent(in) :: A(3, 3)

    determinant_3d_integer = A(1, 1) * A(2, 2) * A(3, 3)  -  A(3, 1) * A(2, 2) * A(1, 3) &
                           + A(1, 2) * A(2, 3) * A(3, 1)  -  A(3, 2) * A(2, 3) * A(1, 1) &
                           + A(1, 3) * A(2, 1) * A(3, 2)  -  A(3, 3) * A(2, 1) * A(1, 2)
  end function determinant_3d_integer

  !> Calculate triple product for three vectors \( \mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3 \in \mathbb{R}^3 \):
  !> \[
  !>   s = \mathbf{a}_{1} \cdot (\mathbf{a}_{2} \times \mathbf{a}_{3}).
  !> \]
  pure function triple_product_vector_inputs(v1, v2, v3) result(t_prod)
    !> Input vectors
    real(dp), intent(in) :: v1(3), v2(3), v3(3)
    real(dp) :: t_prod
    t_prod = determinant_3d_real_dp(reshape([v1, v2, v3], [3, 3]))
  end function

  !> Calculate the cross product between two three dimensional real vectors 
  !> \\(mathbf{b}\) and  \\(mathbf{c}\):
  !> \[ 
  !>   \mathbf{a} = \mathbf{b} \times \mathbf{c}
  !>              = \epsilon_{ijk}b_j c_k
  !> \]
  pure function cross_product(b, c) result(a)
    !> input vectors
    real(dp), intent(in) :: b(3), c(3)
    
    real(dp) :: a(3)

    a(1) = b(2) * c(3) - b(3) * c(2)
    a(2) = b(3) * c(1) - b(1) * c(3)
    a(3)=  b(1) * c(2) - b(2) * c(1)
  end function cross_product

  !> Invert real \(3 \times 3 \) matrix \(A\) such that
  !> \[
  !>   A^{-1} \cdot A = A \cdot A^{-1} = I_{3 \times 3},
  !> \]
  !> where \(I_{3 \times 3}\) is the identity matrix.
  pure function inverse_3d(A) result(inverse_A)
    !> input matrix
    real(dp), intent(in) :: A(3, 3)
    
    real(dp) :: inverse_A(3, 3)

    real(dp) :: determinant

    determinant = determinant_3d(A)

    if (all_zero(determinant)) error stop 'The determinant is too close to zero.'

    inverse_A(1, 1) = (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) / determinant
    inverse_A(1, 2) = (A(1, 3) * A(3, 2) - A(1, 2) * A(3, 3)) / determinant
    inverse_A(1, 3) = (A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)) / determinant

    inverse_A(2, 1) = (A(2, 3) * A(3, 1) - A(2, 1) * A(3, 3)) / determinant
    inverse_A(2, 2) = (A(1, 1) * A(3, 3) - A(1, 3) * A(3, 1)) / determinant
    inverse_A(2, 3) = (A(1, 3) * A(2, 1) - A(1, 1) * A(2, 3)) / determinant
    
    inverse_A(3, 1) = (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1)) / determinant
    inverse_A(3, 2) = (A(1, 2) * A(3, 1) - A(1, 1) * A(3, 2)) / determinant
    inverse_A(3, 3) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) / determinant
  end function inverse_3d

end module linear_algebra_3d