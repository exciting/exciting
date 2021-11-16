module linear_algebra_3d
  use asserts, only: assert
  use precision, only: dp
  use constants, only: pi

  private
  public :: determinant_3d, cross_product, inverse_3d, triple_product, volume_parallelepiped, reciprocal_lattice

  !> Volume tolerance for [[reciprocal_lattice]]. If a parallelpiped is smaller than this,
  !> the volume is taken as zero.
  real(dp), parameter :: v_tol = 1.e-8_dp

  contains 

  ! Calculate the determinant of a matrix \( \mathbf{A} \).
  real(dp) function determinant_3d(A)
    !> Input matrix
    real(dp) :: A(3, 3)

    determinant_3d = A(1, 1) * A(2, 2) * A(3, 3)  -  A(3, 1) * A(2, 2) * A(1, 3) &
                   + A(1, 2) * A(2, 3) * A(3, 1)  -  A(3, 2) * A(2, 3) * A(1, 1) &
                   + A(1, 3) * A(2, 1) * A(3, 2)  -  A(3, 3) * A(2, 1) * A(1, 2)
  end function determinant_3d


  !> Calculate the cross product between two three dimensional real vectors 
  !> \\(mathbf{b}\) and  \\(mathbf{c}\):
  !> \[ 
  !>   \mathbf{a} = \mathbf{b} \times \mathbf{c}
  !>              = \epsilon_{ijk}b_j c_k
  !> \]
  function cross_product(b, c) result(a)

    !> input vectors
    real(dp), intent(in) :: b(3), c(3)
    !> cross product of b and c
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
  function inverse_3d(A) result(inverse_A)
    !> input matrix
    real(dp), intent(in) :: A(3, 3)
    !> inverse of input matrix
    real(dp) :: inverse_A(3, 3)

    real(dp) :: determinant

    determinant = determinant_3d(A)

    call assert(determinant /= 0, 'det A is zero.')

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


  !> Calculate triple product for lattice vectors expected column-wise, 
  !> with the column vectors \( \mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3 \in \mathbb{R}^3}, defined as
  !> \[ 
  !>   s = \mathbf{a}_{1} \cdot (\mathbf{a}_{2} \times \mathbf{a}_{3}),
  !> \]
  !> where
  !> \[ 
  !>   \mathbf{a}_{1} , \mathbf{a}_{2} , \mathbf{a}_{3}
  !> \]
  !> are the lattice vectors.
  real(dp) function triple_product(lattice_vec)
    !> lattice vectors
    real(dp), intent(in) :: lattice_vec(3, 3)

    real(dp) :: b_cross_c(3)

    b_cross_c = cross_product(lattice_vec(:, 2), lattice_vec(:, 3))
    
    triple_product = dot_product(lattice_vec(:, 1), b_cross_c)
  end function triple_product


  !> Calculate the volume of a parallepid, spanned by lattice vectors expected column-wise 
  !> as absolute value of the result of [[triple_product]].
  real(dp) function volume_parallelepiped(lattice_vec)
    !> lattice vectors
    real(dp), intent(in) :: lattice_vec(3, 3)

    volume_parallelepiped = abs(triple_product(lattice_vec))
  end function volume_parallelepiped


  !> Calculate reciprocal lattice vectors from real-space lattice vectors expected column-wise,
  !> \[ 
  !>   \left[\mathbf{b}_{1} \mathbf{b}_{2} \mathbf{b}_{3}\right]^{\top}=
  !>     2 \pi\left[\mathbf{a}_{1} \mathbf{a}_{2} \mathbf{a}_{3}\right]^{-1}
  !> \]
  !> \]
  !> are the real-space lattice vectors and
  !> \[ 
  !>   \mathbf{b}_{1} , \mathbf{b}_{2} , \mathbf{b}_{3} 
  !> \]
  !> are the reciprocal lattice vectors returned column-wise.
  function reciprocal_lattice(lattice) result(rec_lattice)
    !> Lattice vectors
    real(dp), intent(in) :: lattice(3, 3)
    
    real(dp) :: rec_lattice(3, 3)

    real(dp) :: omega

    omega = triple_product(lattice)

    call assert(abs(triple_product(lattice)) > v_tol, 'Unit cell volume is too small.')

    rec_lattice(:, 1) = 2 * pi / omega * cross_product(lattice(:, 2), lattice(:, 3))
    rec_lattice(:, 2) = 2 * pi / omega * cross_product(lattice(:, 3), lattice(:, 1))
    rec_lattice(:, 3) = 2 * pi / omega * cross_product(lattice(:, 1), lattice(:, 2))
  end function reciprocal_lattice

end module linear_algebra_3d