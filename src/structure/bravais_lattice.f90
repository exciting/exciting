!> Contains functions for generating 3d Bravais lattices.
!> The implementation and documentation is based on the definitions in the 
!> [AFLOW enceclopedia](http://aflowlib.org/prototype-encyclopedia/).
module bravais_lattice
  use precision, only: dp
  use constants, only: pi
  use asserts, only: assert

  implicit none

  private
  public :: simple_cubic, &
            body_centered_cubic, &
            face_centered_cubic, &

            hexagonal, &
            graphene, &
            rhombohedral_hex_setting, &
            rhombohedral_rhom_setting, &

            simple_tetragonal, &
            body_centered_tetragonal, &

            simple_orthorhombic, &
            base_centered_orthorhombic_A, &
            base_centered_orthorhombic_C, &
            body_centered_orthorhombic, &
            face_centered_orthorhombic, &

            simple_monoclinic, &
            base_centered_monoclinic_unquie_axis_c, &
              
            triclinic

  real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
                          
contains

  ! Cubic

  !> Simple cubic conventional lattice vectors.
  !>
  !> The cubic crystal system is defined as having the symmetry of a cube: the conventional unit 
  !> cell can be rotated by \(\pi / 2\) about any axis, or by \(\pi\) around an axis running through the center 
  !> of two opposing cube edges, or by \( 2 \pi / 3 \) around a body diagonal, and retain the same shape.
  !>
  !> This definition and documentation follows 
  !> [AFLOW cubic lattice](http://aflowlib.org/prototype-encyclopedia/cubic_lattice.html).
  function simple_cubic(a) result(lattice)
    !> lattice constant
    real(dp), intent(in) :: a
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)

    lattice = a * transpose(reshape([1, 0, 0, &
                                     0, 1, 0, &
                                     0, 0, 1], [3, 3]))
  end function simple_cubic


  !> Face centered cubic lattice.
  !>
  !> The face-centered cubic lattice has the same periodicity as its simple cubic parent with the 
  !> addition of a translation from one corner of the cube to the center of any face.
  !>
  !> This definition and documentation follows 
  !> [AFLOW cubic lattice](http://aflowlib.org/prototype-encyclopedia/cubic_lattice.html).
  function face_centered_cubic(a) result(lattice)
    !> lattice constant
    real(dp), intent(in) :: a
    
    real(dp) :: lattice(3, 3)

    lattice = 0.5_dp * a * transpose(reshape([0, 1, 1, &
                                              1, 0, 1, &
                                              1, 1, 0], [3, 3]))
  end function


  !> Body-centered cubic lattice
  !>
  !> The body-centered cubic crystal has the same periodicity as its parent with the addition 
  !> of a translation from one corner of the cube to its center.
  !>
  !> This definition and documentation follows 
  !> [AFLOW cubic lattice](http://aflowlib.org/prototype-encyclopedia/cubic_lattice.html).
  function body_centered_cubic(a) result(lattice)
    !> lattice constant
    real(dp), intent(in) :: a
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)

    lattice = 0.5_dp * a * transpose(reshape([-1, 1, 1, &
                                               1,-1, 1, &
                                               1, 1,-1], [3, 3]))
  end function body_centered_cubic

  
  ! Trigonal/Hexagonal

  !> Hexagonal conventional lattice vectors
  !>
  !> Also referred to as the simple trigonal lattice.
  !> Shares conventional vectors but not point operations with the hexagonal crystal system.
  !> Note, one could equivalently define \(\mathbf{a_1} = {1,0,0}\)
  !>
  !> This definition and documentation follows 
  !> [AFLOW trigonal crystal system](http://aflowlib.org/prototype-encyclopedia/trigonal_lattice.html).
  function hexagonal(a, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(c > 0)

    lattice = a * transpose(reshape([        0.5_dp,         0.5_dp, 0.0_dp, &
                                    -0.5_dp * sqrt3, 0.5_dp * sqrt3, 0.0_dp, &
                                             0.0_dp,         0.0_dp,    c/a], [3, 3]))
  end function


  !> Graphene lattice 
  function graphene(a, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(c > 0)

    lattice = a * transpose(reshape([1.0_dp,         0.5_dp, 0.0_dp, &
                                     0.0_dp, 0.5_dp * sqrt3, 0.0_dp, &
                                     0.0_dp,         0.0_dp,    c/a], [3, 3]))
  end function


  !> Rhombohedral lattice vectors (hexagonal setting)
  !>
  !> One can define the rhombohedral lattice in two ways: as a hexagonal lattice
  !> with additional translational vectors (hexagonal setting), or as a simple
  !> lattice with conventional vectors of equal length making equal angles with one
  !> another (rhombohedral setting).
  !>
  !> This definition and documentation follows 
  !> [AFLOW trigonal crystal system](http://aflowlib.org/prototype-encyclopedia/trigonal_lattice.html).
  function rhombohedral_hex_setting(a, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, c
    
    real(dp) :: lattice(3, 3)

    real(dp) :: x, y, z 

    call assert(a > 0)
    call assert(c > 0)

    x = 0.5_dp * a
    y = 0.5_dp * a / sqrt3
    z = c / 3.0_dp

    lattice = transpose(reshape([ x,     0.0_dp, -x, &
                                 -y, 2.0_dp * y, -y, &
                                  z,          z,  z], [3, 3]))
  end function


  !> Rhombohedral lattice vectors (rhom setting)
  !>
  !> One can define the rhombohedral lattice in two ways: as a hexagonal lattice
  !> with additional translational vectors (hexagonal setting), or as a simple
  !> lattice with conventional vectors of equal length making equal angles with one
  !> another (rhombohedral setting).
  !>
  !> This definition and documentation follows 
  !> [AFLOW trigonal crystal system](http://aflowlib.org/prototype-encyclopedia/trigonal_lattice.html).
  !>
  !> Angles are given in radian.
  function rhombohedral_rhom_setting(a, alpha) result(lattice)
    !> lattice constant
    real(dp), intent(in) :: a
    !> lattice angle in degree
    real(dp), intent(in) :: alpha
    
    real(dp) :: lattice(3, 3)

    real(dp) :: x, y, cx, cz

    call assert(a > 0)
    call assert(alpha > 0)

    x = sin(0.5_dp * alpha)
    y = sin(0.5_dp * alpha) / sqrt3
    cz = sqrt((4._dp * cos(0.5 * alpha)**2 - 1._dp)) / sqrt3;

    lattice = a * transpose(reshape([ x,     0.0_dp, -x, &
                                     -y, 2.0_dp * y, -y, &
                                     cz,         cz, cz], [3, 3]))
  end function


  ! Tetragonal

  !> Simple tetragonal conventional lattice vectors.
  !> 
  !> In the tetragonal lattice system, the conventional unit cell is a parallelepiped
  !> but two sides are equal such that 'a = b' and \(c \neq a\) while \(\alpha = \beta = \gamma = \pi / 2 \),
  !> and this is a special case of the orthorhombic system.
  !>
  !> This definition and documentation follows 
  !> [AFLOW tetragonal crystal system](http://aflowlib.org/prototype-encyclopedia/tetragonal_lattice.html).
  function simple_tetragonal(a, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(c > 0)

    lattice = transpose(reshape([     a, 0.0_dp, 0.0_dp, &
                                 0.0_dp,      a, 0.0_dp, &
                                 0.0_dp, 0.0_dp,      c], [3, 3]))
  end function


  !> Body-centered tetragonal lattice vectors.
  !>
  !> The body-centered tetragonal lattice system has the same point group and translational
  !> symmetry as  the simple tetragonal system, with addition of a translation to the center
  !> of the parallelpiped.
  !>
  !> This definition and documentation follows 
  !> [AFLOW tetragonal crystal system](http://aflowlib.org/prototype-encyclopedia/tetragonal_lattice.html).
  function body_centered_tetragonal(a, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(c > 0)

    lattice = 0.5_dp * transpose(reshape([-a,  a,  a, &
                                           a, -a,  a, &
                                           c,  c, -c], [3, 3]))
  end function


  ! Orthorhombic

  !> Simple orthorhombic conventional lattoce vectors. 
  
  !> In the orthorhombic system the conventional unit cell is a parallelpiped defined by three 
  !> mutionally orthogonal vectors of unequal length so that 
  !> \( a \neq b \neq c \), but \( \alpha = \beta = \gamma \). It is a limiting case of the 
  !> monoclinic crystal with \( \beta \rightarrow \pi / 2 \).
  !>
  !> This definition and documentation follows 
  !> [AFLOW orthorhombic crystal system](http://aflowlib.org/prototype-encyclopedia/orthorhombic_lattice.html).
  !>
  !> Constant 'a' should be smaller than 'b' and 'c'.
  function simple_orthorhombic(a, b, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)

    call assert(a < b)
    call assert(a < c)

    lattice = transpose(reshape([     a, 0.0_dp, 0.0_dp, &
                                 0.0_dp,      b, 0.0_dp, &
                                 0.0_dp, 0.0_dp,      c], [3, 3]))
  end function


  !> Base-centered orthorhombic lattice vectors in the C setting
  !>
  !> Base-centered cell allows for a translation in one of the basal planes.
  !> Space groups ending with C put the translation in the 'a-b' plane.
  !>
  !> This definition and documentation follows 
  !> [AFLOW orthorhombic crystal system](http://aflowlib.org/prototype-encyclopedia/orthorhombic_lattice.html).
  function base_centered_orthorhombic_C(a, b, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)

    lattice = 0.5_dp * transpose(reshape([     a,      a,     0.0_dp, &
                                              -b,      b,     0.0_dp, &
                                          0.0_dp, 0.0_dp, 2.0_dp * c], [3, 3]))
  end function


  !> Base-centered orthorhombic lattice vectors in the A setting
  !>
  !> base-centered cell allows for a translation in one of the basal planes.
  !> Space groups ending with A put the translation in the 'b−c' plane.
  !>
  !> This definition and documentation follows 
  !> [AFLOW orthorhombic crystal system](http://aflowlib.org/prototype-encyclopedia/orthorhombic_lattice.html).
  !> by rotating the system about \(\pi / 2\) around \(\mathbf{b}\).
  function base_centered_orthorhombic_A(a, b, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)

    lattice = 0.5 * transpose(reshape([2.0_dp * a,  0.0_dp,  0.0_dp, &
                                           0.0_dp,       b,       b, &
                                           0.0_dp,      -c,       c], [3, 3]))
  end function


  !> Body-centered orthorhombic
  !>
  !> The body-centered orthorhombic lattice has the same point group and translational 
  !> symmetry as the simple orthorhombic system, with the addition of a translation 
  !> to the center of the parallelepiped.
  !>
  !> This definition and documentation follows 
  !> [AFLOW orthorhombic crystal system](http://aflowlib.org/prototype-encyclopedia/orthorhombic_lattice.html).
  function body_centered_orthorhombic(a, b, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)

    lattice = 0.5 * transpose(reshape([-a,  a,  a, &
                                        b, -b,  b, &
                                        c,  c, -c], [3, 3]))
  end function


  !> Face-centered orthorhombic
  !>
  !> While the base-centered monoclinic lattice allows translations to one base plane, 
  !> the face-centered orthorhombic lattice allows translations to any of the base planes. 
  !> Our standard choice for the primitive vectors of this system are given by 
  !>
  !> This definition and documentation follows 
  !> [AFLOW orthorhombic crystal system](http://aflowlib.org/prototype-encyclopedia/orthorhombic_lattice.html).
  function face_centered_orthorhombic(a, b, c) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    
    real(dp) :: lattice(3, 3)

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)

    lattice = 0.5 * transpose(reshape([0.0_dp,      a,      a, &
                                            b, 0.0_dp,      b, &
                                            c,      c, 0.0_dp], [3, 3]))
  end function


  ! Monoclinic

  !> Simple monoclinic conventional lattice vectors
  !>
  !> In the monoclinic crystal system, the conventional unit cell is defined by
  !> primitive vectors of arbitrary length, where one of the vectors is perpendicular to
  !> the other two. Modern convention chooses this vector to be \( \mathbf{b} \) such that
  !> \( a \neq b \neq c), and \( \alpha = \gamma = \pi / 2, \beta \neq \pi / 2 \).
  !> 
  !> This definition and documentation follows 
  !> [AFLOW monoclinic crystal system](http://aflowlib.org/prototype-encyclopedia/monoclinic_lattice.html).
  !>
  !> Angles are given in radian.
  function simple_monoclinic(a, b, c, beta) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    !> lattice angle in rad
    real(dp), intent(in) :: beta
    
    real(dp) :: lattice(3, 3)

    real(dp) :: cx, cz

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)
    call assert(beta > 0)

    cx = c * cos(beta)
    cz = c * sin(beta)

    lattice = transpose(reshape([     a, 0.0_dp,     cx, &
                                 0.0_dp,      b, 0.0_dp, &
                                 0.0_dp, 0.0_dp,     cz], [3, 3]))
  end function


  !> Base-centered monoclinic lattice vectors
  !>
  !> The base-centered monoclinic lattice is in the same crystal system as the monoclinic lattice, 
  !> but its periodicity allows an additional translation in the plane defined by \( \mathbf{a1} \) 
  !> and \( \mathbf{a2} \). 
  !> 
  !> This definition and documentation follows 
  !> [AFLOW monoclinic crystal system](http://aflowlib.org/prototype-encyclopedia/monoclinic_lattice.html).
  !>
  !> Angles are given in radian.
  function base_centered_monoclinic_unquie_axis_c(a, b, c, beta) result(lattice)
    !> lattice constants
    real(dp), intent(in) :: a, b, c
    !> lattice angle in rad
    real(dp), intent(in) :: beta
    
    real(dp) :: lattice(3, 3)

    real(dp) :: cx, cz

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)
    call assert(beta > 0)
    
    cx = c * cos(beta)
    cz = c * sin(beta)

    lattice = transpose(reshape([0.5_dp * a, 0.5_dp * a,     cx, &
                                -0.5_dp * b, 0.5_dp * b, 0.0_dp, &
                                 0.0_dp,     0.0_dp,         cz], [3, 3]))
  end function


  ! Triclinic

  !> Triclinic conventional lattice vectors
  !>
  !> Triclinic is the most general crystal system. All other crystal systems can be considered special cases of the triclinic. 
  !> The primitive vectors are also completely general: their lenghts \( a, b, c \) and angles 
  !> \(\alpha, \beta, \gamma) may have arbitrary values. The triclinic system has one Bravais lattice, 
  !> which is also the conventional lattice for this system. 
  !>
  !> This definition and documentation follows 
  !> [AFLOW triclinic lattice](http://aflowlib.org/prototype-encyclopedia/triclinic_lattice.html).
  !>
  !> Angles are given in radian.
  function triclinic(a, b, c, alpha, beta, gamma) result(lattice)
    !> Lattice constants
    real(dp), intent(in) :: a, b, c
    !> Lattice angles in rad
    real(dp), intent(in) :: alpha, beta, gamma
    
    real(dp) :: lattice(3, 3)

    real(dp) :: bx, by, cx, cy, cz

    call assert(a > 0)
    call assert(b > 0)
    call assert(c > 0)
    call assert(alpha > 0); call assert(alpha < pi / 2)
    call assert(beta > 0); call assert(beta < pi / 2)
    call assert(gamma > 0); call assert(gamma < pi / 2)

    bx = b * cos(gamma)
    by = b * sin(gamma)
    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)

    call assert(c * c > cx * cx + cy * cy)
    cz = sqrt(c*c - cx*cx - cy*cy)
    
    lattice = transpose(reshape([     a,     bx, cx, &
                                 0.0_dp,     by, cy, &
                                 0.0_dp, 0.0_dp, cz], [3, 3]))
  end function

end module bravais_lattice