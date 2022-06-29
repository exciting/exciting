module linear_algebra_3d_test
  use precision, only: dp
  use constants, only: pi
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use math_utils, only: identity_real_dp, all_close, all_zero
  use bravais_lattice, only: simple_cubic, face_centered_cubic, body_centered_cubic
  use linear_algebra_3d, only: determinant_3d, &
                               cross_product, &
                               inverse_3d, &
                               triple_product

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
    integer, parameter :: n_assertions = 23

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    
    call test_determinant_3d(test_report)

    call test_cross_product(test_report)

    call test_inverse_3d(test_report)

    call test_triple_product(test_report)

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

    real(dp) :: test_vectors(3, 3)

    ! v_1 * (v_2 x v_3)

    test_vectors = transpose(reshape([1._dp, 0._dp, 0._dp, &
                                      0._dp, 1._dp, 0._dp, &
                                      0._dp, 0._dp, 1._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  1._dp), &
                            'Triple product is not 1.0 for 3d identity as input.')

    test_vectors = transpose(reshape([-1._dp, 0._dp, 0._dp, &
                                       0._dp, 1._dp, 0._dp, &
                                       0._dp, 0._dp, 1._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  -1._dp), &
                            'Triple product is not -1.0 for orthogonal unit vecotors as input, where the first vector &
                            is negative.')

    test_vectors = transpose(reshape([1._dp,  0._dp, 0._dp, &
                                      0._dp, -1._dp, 0._dp, &
                                      0._dp,  0._dp, 1._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  -1._dp), &
                            'Triple product is not -1.0 for orthogonal unit vecotors as input, where the second vector &
                            is negative.')

    test_vectors = transpose(reshape([1._dp,  0._dp,  0._dp, &
                                      0._dp,  1._dp,  0._dp, &
                                      0._dp,  0._dp, -1._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  -1._dp), &
                            'Triple product is not -1.0 for orthogonal unit vecotors as input, where the third vector &
                            is negative.')

    test_vectors = transpose(reshape([-1._dp,  0._dp,  0._dp, &
                                       0._dp,  1._dp,  0._dp, &
                                       0._dp,  0._dp, -1._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  1._dp), &
                            'Triple product is not 1.0 for orthogonal unit vecotors as input, where the first and the &
                             third vectors are negative.')

    test_vectors = transpose(reshape([1._dp,  0._dp, 0._dp, &
                                      0._dp,  1._dp, 1._dp, &
                                      0._dp,  0._dp, 0._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  0._dp), &
                            'Triple product is not 0.0 where the second and the third vectors are linearly &
                            dependendent.')

    test_vectors = transpose(reshape([0._dp,  0._dp, 0._dp, &
                                      1._dp,  1._dp, 0._dp, &
                                      0._dp,  0._dp, 1._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  0._dp), &
                            'Triple product is not 0.0 where the first and the second vectors are linearly &
                            dependendent.')

    test_vectors = transpose(reshape([1._dp,  0._dp, 1._dp, &
                                      0._dp,  1._dp, 0._dp, &
                                      0._dp,  0._dp, 0._dp ], [3, 3]))
    call test_report%assert(all_close(triple_product(test_vectors),  0._dp), &
                            'Triple product is not 0.0 where the first and the third vectors are linearly &
                            dependendent.')

    test_vectors = transpose(reshape([ 1._dp,  -0.23_dp, 1.43_dp, &
                                       5._dp,    1.1_dp, -2.3_dp, &
                                      10._dp,    12._dp,   3._dp ], [3, 3]))

    call test_report%assert(all_close(triple_product(test_vectors), &
                                      triple_product(test_vectors(:, 1), test_vectors(:, 2), test_vectors(:, 3))), &
                            'Triple product interfaces for matrix input and for vector input give not the same result.')
  end subroutine test_triple_product

end module linear_algebra_3d_test