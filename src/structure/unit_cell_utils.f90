module unit_cell_utils
  use asserts, only: assert
  use precision, only: dp
  use constants, only: pi
  use math_utils, only: all_zero
  use linear_algebra_3d, only: triple_product, inverse_3d

  implicit none

  private
  public :: volume_parallelepiped, &
            reciprocal_lattice

  contains

  !> Calculate the volume of a parallelepiped, spanned by lattice vectors expected column-wise
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
  !> are the real-space lattice vectors and
  !> \[ 
  !>   \mathbf{b}_{1} , \mathbf{b}_{2} , \mathbf{b}_{3} 
  !> \]
  !> are the reciprocal lattice vectors returned column-wise.
  pure function reciprocal_lattice(lattice) result(rec_lattice)
    !> Lattice vectors
    real(dp), intent(in) :: lattice(3, 3)
    
    real(dp) :: rec_lattice(3, 3)

    real(dp) :: omega

    omega = triple_product(lattice)

    if(all_zero(omega)) error stop 'The unit cell volume is too close to zero.'

    rec_lattice = 2 * pi * transpose(inverse_3d(lattice))
  end function reciprocal_lattice




end module unit_cell_utils