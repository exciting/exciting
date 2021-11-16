module linear_algebra_3d_test
  use precision, only: dp
  use constants, only: pi
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use math_utils, only: identity_real_dp, &
                        all_close, &
                        all_zero
  use linear_algebra_3d, only: determinant_3d, &
                               cross_product, &
                               inverse_3d, &
                               triple_product, &
                               volume_parallelepiped, &
                               reciprocal_lattice

  use bravais_lattice, only: simple_cubic, &
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

  implicit none

  private
  public :: linear_algebra_3d_test_driver

  real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
    
  
contains

  !> Run tests for math tools
  subroutine linear_algebra_3d_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure  
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 48

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    
    call test_determinant_3d(test_report)

    call test_cross_product(test_report)

    call test_inverse_3d(test_report)

    call test_triple_product(test_report)

    call test_reciprocal_lattice(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('linear_algebra_3d', kill_on_failure)
    else
      call test_report%report('linear_algebra_3d')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine linear_algebra_3d_test_driver


  !> Test determinant_3d
  subroutine test_determinant_3d(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    !> test input
    real(dp) :: test_matrix(3, 3) = reshape([ 1.0_dp, 0.0_dp, 6.0_dp, &
                                              0.0_dp, 3.0_dp, 0.0_dp, &
                                              5.0_dp, 0.0_dp, 4.0_dp ], [ 3, 3 ])
           
    call test_report%assert(all_close(determinant_3d(test_matrix), -78.0_dp), &
               'Test determinant_3d for real bijective matrix (the deteminant is not zero). &
               Expected: -78.0')

    test_matrix = reshape([ 1.0, 2.0, 0.0, &
                            3.0, 4.0, 0.0, &
                            0.0, 0.0, 0.0 ], [ 3, 3 ])
    call test_report%assert(all_close(determinant_3d(test_matrix), 0.0_dp), &
                            'Test determinant_3d for real non-bijective matrix (the determinant_3d is zero). &
                            Expected: 0.0')
  end subroutine test_determinant_3d


  !> Test cross_product
  subroutine test_cross_product(test_report)
    !> test object
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: a(3), b(3), c(3), zero(3), &
                a_cross_b(3), b_cross_a(3), b_cross_c(3), c_cross_a(3), & 
                norm_a, norm_b, norm_a_cross_b

    call test_report%assert(all_close(cross_product([3.0_dp, 0.0_dp, 0.0_dp], &
                                                    [0.0_dp, 2.0_dp, 0.0_dp]), &
                                                    [0.0_dp, 0.0_dp, 6.0_dp]), &
                            'Test cross_product for unit vector in x and y direction. &
                            Expected: Vector in z direction.')

    call test_report%assert(all_close(cross_product([3.0_dp, 0.0_dp, 0.0_dp], &
                                                    [2.0_dp, 0.0_dp, 0.0_dp]), &
                                                    [0.0_dp, 0.0_dp, 0.0_dp]), &
                            'Test cross_product for collinear vectors. &
                            Expected: Zero vector.')

    ! set up variables for testing
    zero = 0.0_dp

    a = [0.23_dp, 22.41_dp, -24.34_dp]
    b = [-21.04_dp, 0.30_dp, 24.12_dp]
    c = [21.12_dp, 0.23_dp, -0.92_dp]

    a_cross_b = cross_product(a, b)
    b_cross_a = cross_product(b, a)
    b_cross_c = cross_product(b, c)
    c_cross_a = cross_product(c, a)

    norm_a_cross_b = norm2(a_cross_b)
    norm_a = norm2(a)
    norm_b = norm2(b)

    call test_report%assert(all_close(dot_product(a_cross_b, a), 0.0_dp), &
                            'Test that a x b is orthogonal to a. &
                            Expected: (a x b) * a = 0')

   
    call test_report%assert(all_close(dot_product(b, a_cross_b), 0.0_dp), &
                            'Test that a x b is orthogonal to b. &
                            Expected: (a x b) * b = 0')

    call test_report%assert(all_close(norm_a_cross_b, sqrt(abs(norm_a**2 * norm_b**2 - dot_product(a, b)**2))), &
                           'Test that a x b has the correct length. &
                           Expected: |a x b|^2 = |a|^2 * |b|^2 - (a * b)^2 ')


    call test_report%assert(all_close(a_cross_b, -1.0_dp * b_cross_a), &
                            'Test that cross_product is anti commutative. &
                            Expected: a x b == - b x a')

    call test_report%assert(all_zero(cross_product(a, b_cross_c) & 
                                   + cross_product(b, c_cross_a) &
                                   + cross_product(c, a_cross_b), 1e-8_dp), &
                           'Test the Jacobi identity for three vectors a, b and c. &
                           Expected: a x (b x c) + b x (c x a) + c x (a x b) = 0_vec')
  end subroutine test_cross_product

  !> Test invert
  subroutine test_inverse_3d(test_report)
    !> test object
    type(unit_test_type), intent(inout) :: test_report
    !> reference matrices
    real(dp) :: sc(3, 3), fcc(3, 3), bcc(3, 3)

    sc = simple_cubic(4.0_dp)

    fcc = face_centered_cubic(3.5_dp)

    bcc = body_centered_cubic(7.5_dp)

    call test_report%assert(all_close(matmul(inverse_3d(sc), sc), identity_real_dp(3)), &
                           'Test inverse_3d for sc lattice vectors by multiplying  &
                           the result of inverse_3d to the matrix itself. &
                           Expected: identity matrix')


    call test_report%assert(all_close(matmul(inverse_3d(fcc), fcc), identity_real_dp(3)), &
                           'Test inverse_3d for fcc lattice by multiplying  &
                           the result of inverse_3d to the matrix itself. &
                           Expected: identity matrix')


    call test_report%assert(all_close(matmul(inverse_3d(bcc), bcc), identity_real_dp(3)), &
                           'Test inverse_3d for bcc lattice by multiplying  &
                           the result of inverse_3d to the matrix itself. &
                           Expected: identity matrix')

  end subroutine test_inverse_3d


  !> Test triple product
  subroutine test_triple_product(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: omega, a, aprime, b, c, alpha, beta, gamma, cx, cy, cz

    real(dp), parameter :: sqrt3 = sqrt(3.0_dp)


    call test_report%assert(all_close(triple_product( 0.5 * body_centered_cubic(1.0_dp)),  6.25e-2_dp), &
                            'Test if triple_product can be positive.')

    call test_report%assert(all_close(triple_product(-0.5 * body_centered_cubic(1.0_dp)),  -6.25e-2_dp), &
                            'Test if triple_product can be negative.')

    ! cubic lattice vectors
    a = 2.0_dp
    omega = volume_parallelepiped(simple_cubic(a))
    call test_report%assert(all_close(omega, a**3), &
                           'Test triple product for simple cubic lattice vectors. &
                           Expected: omega = a**3')

    omega = volume_parallelepiped(face_centered_cubic(a))
    call test_report%assert(all_close(omega, a**3 / 4), &
                           'Test triple product for face-centered cubic lattice vectors. &
                           Expected: a**3 / 4')

    omega = volume_parallelepiped(body_centered_cubic(a))
    call test_report%assert(all_close(omega, a**3 / 2), &
                           'Test triple product for body-centered cubic lattice vectors. &
                           Expected: a**3 / 2')

    ! Trigonal lattice vectors
    a = 2.5; c = 1.4
    omega = volume_parallelepiped(hexagonal(a, c))
    call test_report%assert(all_close(omega, sqrt3 / 2 * a**2 * c), &
                           'Test triple product for hexagonal lattice vectors. &
                           Expected: omega = sqrt3 / 2 * a**2 * c')

    omega = volume_parallelepiped(rhombohedral_hex_setting(a, c))
    call test_report%assert(all_close(omega, a**2 * c / (2 * sqrt3)), &
                           'Test triple product for rhombohedral lattice vectors (hexagonal setting). &
                           Expected: omega = a**2 * c / (2 * sqrt3)')

    ! The formulas for the rhombohedral setting are from AFLOW (see doc)
    aprime = sqrt(a**2 / 3 + c**2 / 9)
    alpha = dacos((2 * c**2 - 3 * a**2) / (2 * (c**2 + 3 * a**2)))
    omega = volume_parallelepiped(rhombohedral_rhom_setting(aprime, alpha))
    call test_report%assert(all_close(omega, a**2 * c / (2 * sqrt3)), &
                           'Test triple product for rhombohedral lattice vectors (rhombohedral setting). &
                           Expected: omega = a**2 * c / (2 * sqrt3)')

    ! Tetragonal
    a = 2.5; c = 1.4
    omega = volume_parallelepiped(simple_tetragonal(a, c))
    call test_report%assert(all_close(omega, a**2 * c), &
                           'Test triple product for simple tetragonal lattice vectors. &
                           Expected: omega = a**2 * c')
                           
    omega = volume_parallelepiped(body_centered_tetragonal(a, c))
    call test_report%assert(all_close(omega, a**2 * c / 2), &
                           'Test triple product for body-centered tetragonal lattice vectors. &
                           Expected: omega = a**2 * c / 2')

    ! Orthorhombic
    a = 1.5; b = 3.1; c = 1.6
    omega = volume_parallelepiped(simple_orthorhombic(a, b, c))
    call test_report%assert(all_close(omega, a * b * c), &
                           'Test triple product for simple orthorhombic lattice vectors. &
                           Expected: omega = a * b * c')
                           
    omega = volume_parallelepiped(base_centered_orthorhombic_C(a, b, c))
    call test_report%assert(all_close(omega, a * b * c / 2), &
                           'Test triple product for base-centered orthorhombic lattice vectors (C setting). &
                           Expected: omega = a * b * c / 2')

    omega = volume_parallelepiped(base_centered_orthorhombic_A(a, b, c))
    call test_report%assert(all_close(omega, a * b * c / 2), &
                           'Test triple product for base-centered orthorhombic lattice vectors (A setting). &
                           Expected: omega = a * b * c / 2')

    omega = volume_parallelepiped(body_centered_orthorhombic(a, b, c))
    call test_report%assert(all_close(omega, a * b * c / 2), &
                           'Test triple product for body-centered orthorhombic lattice vectors. &
                           Expected: omega = a * b * c / 2')

    omega = volume_parallelepiped(face_centered_orthorhombic(a, b, c))
    call test_report%assert(all_close(omega, a * b * c / 4), &
                           'Test triple product for face-centered orthorhombic lattice vectors. &
                           Expected: omega = a * b * c / 4')

    ! Monoclinic
    a = 1.5; b = 3.1; c = 1.4
    beta = 1.34 ! angles in radian
    omega = volume_parallelepiped(simple_monoclinic(a, b, c, beta))
    call test_report%assert(all_close(omega, a * b * c * sin(beta)), &
                           'Test triple product for simple monoclinic lattice vectors. &
                           Expected: omega = a * b * c * sin(beta)')

    omega = volume_parallelepiped(base_centered_monoclinic_unquie_axis_c(a, b, c, beta))
    call test_report%assert(all_close(omega, a * b * c * sin(beta) / 2), &
                           'Test triple product for base-centered monoclinic lattice vectors. &
                           Expected: omega = a * b * c * sin(beta) / 2')

    ! Triclinic
    a = 1.5; b = 3.1; c = 1.6 
    alpha = 0.92_dp; beta = 1.34_dp; gamma = 1.23_dp 

    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)
    cz = sqrt(c**2 - cx**2 - cy**2)

    omega = volume_parallelepiped(triclinic(a, b, c, alpha, beta, gamma))
    call test_report%assert(all_close(omega, a * b * cz * sin(gamma)), &
                           'Test triple product for simple monoclinic lattice vectors. &
                           Expected: omega = a * b * cz * sin(gamma)')

    
  end subroutine test_triple_product

  !> Test reciprocal lattice          
  subroutine test_reciprocal_lattice(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report
    !> inverse real spave lattice vectors
    real(dp) :: two_pi_id(3, 3), lattice(3, 3), rec_lattice(3, 3)
    real(dp) :: omega, a, aprime, b, c, alpha, beta, gamma

    
    two_pi_id = 2 * pi * identity_real_dp(3)
    
    ! Cubic lattice vectors
    a = 254.2345_dp
    lattice = simple_cubic(a)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for simple cubic lattice. &
                            Expected: lattice**T * rec_lattic = 2 * pi * delta.')
    
    call test_report%assert(all_close(rec_lattice, simple_cubic(2.0_dp * pi / a )), &
                            'Test reciprocal_lattice for simple cubic lattice. &
                            Expected: simple cubic lattice with scaled lattice constant.')

    lattice = face_centered_cubic(a)
    omega = abs(triple_product(lattice))
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for face-centered cubic lattice. &
                            Expected: lattice**T * rec_lattic = 2 * pi * delta.')

    call test_report%assert(all_close(rec_lattice, body_centered_cubic(4.0_dp * pi / a )), &
                            'Test reciprocal_lattice for face-centered cubic lattice. &
                            Expected: body centered cubic with scaled lattice constant.')

    lattice = body_centered_cubic(a)
    omega = abs(triple_product(lattice))
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for body-centered cubic lattice. &
                            Expected: lattice**T * rec_lattic = 2 * pi * delta.')

    call test_report%assert(all_close(rec_lattice, face_centered_cubic(4.0_dp * pi / a)), &
                            'Test reciprocal_lattice for body-centered cubic lattice. &
                            Expected: face centered cubic with scaled lattice constant.')

    ! Hexagonal lattice vectors
    a = 21.345_dp; c = 142.24_dp
    lattice = hexagonal(a, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for hexagonal lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')


    lattice = rhombohedral_hex_setting(a, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for hexagonal lattice (hex. setting). &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    aprime = sqrt(a**2 / 3 + c**2 / 9)
    alpha = dacos((2 * c**2 - 3 * a**2) / (2 * (c**2 + 3 * a**2)))
    lattice = rhombohedral_rhom_setting(aprime, alpha)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for hexagonal lattice (rhom. setting). &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    ! Tetragonal
    a = 3.24_dp; c = 4.325_dp
    lattice = simple_tetragonal(a, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for simple tetragonal lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    lattice = body_centered_tetragonal(a, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for body-centered tetragonal lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    ! Orthorhombic
    a = 23.24_dp; b = 34.24145_dp; c = 25.325_dp
    lattice = base_centered_orthorhombic_C(a, b, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for base-centered orthorhombic lattice (C setting). &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    lattice = base_centered_orthorhombic_A(a, b, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for base-centered orthorhombic lattice (A setting). &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    lattice = body_centered_orthorhombic(a, b, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for body-centered orthorhombic lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    lattice = face_centered_orthorhombic(a, b, c)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for face-centered orthorhombic lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    ! Monoclinic
    a = 0.234_dp; b = 2.3412_dp; c = 421._dp 
    beta = 1.243_dp ! angles are in radian
    lattice = simple_monoclinic(a, b, c, beta)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for simple monoclinic lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    lattice = base_centered_monoclinic_unquie_axis_c(a, b, c, beta)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for base-centered monoclinic lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.')

    ! Triclinic
    a = 1.5132_dp; b = 1.1_dp; c = 13.64543_dp
    alpha = 0.32492; beta = 1.144_dp; gamma = 1.123_dp
    lattice = triclinic(a, b, c, alpha, beta, gamma)
    rec_lattice = reciprocal_lattice(lattice)
    call test_report%assert(all_close(matmul(transpose(lattice), rec_lattice), two_pi_id), & 
                            'Test reciprocal_lattice for triclinic lattice. &
                            Expected: lattice**T * rec_lattice = 2 * pi * delta.') 
    
                
  end subroutine test_reciprocal_lattice
end module linear_algebra_3d_test