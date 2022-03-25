!> Module for unit tests for for the functions in math_utils.

module math_utils_test
  use precision, only: dp
  use constants, only: zone, zzero
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type
  use math_utils, only: all_close, &
                        all_zero, &
                        diag, &
                        is_square, &
                        is_hermitian, &
                        is_unitary, &
                        kronecker_product, &
                        determinant, &
                        permanent, &
                        mod1, &
                        shuffle_vector, &
                        mask_vector, &
                        round_down, &
                        distance_matrix, &
                        outer_sum, &
                        get_subinterval_indices 
  use mock_arrays, only: real_symmetric_matrix_5x5, &
                         real_orthogonal_matrix_5x5, &
                         real_matrix_5x7, &
                         real_rank_4_matrix_5x7, &
                         complex_hermitian_matrix_5x5, &
                         complex_unitary_matrix_5x5, &
                         complex_matrix_5x7
                         

  implicit none

  real(dp), parameter :: identity_3d(3, 3) = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                                                      0.0_dp, 1.0_dp, 0.0_dp, &
                                                      0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
  private
  public :: math_utils_test_driver

contains

  !> Run tests for math tools
  subroutine math_utils_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 73

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests

    call test_all_close(test_report)

    call test_all_zero(test_report)

    call test_is_square(test_report)

    call test_is_hermitian(test_report)

    call test_is_unitary(test_report)

    call test_diag(test_report)
    
    call test_kronecker_product(test_report)

    call test_determinant(test_report)

    call test_mod1(test_report)

    call test_mask_vector(test_report)

    call test_round_down(test_report)

    call test_distance_matrix(test_report)

    call test_outer_sum(test_report)

    call test_get_subinterval_indices(test_report)

  
    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('math_utils', kill_on_failure)
    else
      call test_report%report('math_utils')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine math_utils_test_driver

  
  !> Test all_close for real and complex arrays
  subroutine test_all_close(test_report)

    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    !> tolerance
    real(dp), parameter :: tol = 1.e-10_dp
    !> Real test input
    real(dp) :: a(3, 4) 
    !> Complex test input
    complex(dp) :: b(3,4)

    a =  1._dp 
    b = cmplx(1.0, 1.0, dp) 

    call test_report%assert(all_close(a, a, tol), &
                  & 'Expected: Real arrays are equal')
    call test_report%assert(.not. all_close(a, a + 2*tol, tol),&
                  & 'Expected: Real arrays are not equal')

    call test_report%assert(all_close(b, b, tol), &
                  & 'Expected: Complex arrays are equal')
    call test_report%assert(.not. all_close(b, b + 2*tol, tol), &
                  & 'Expected: Complex arrays are not equal')

  end subroutine test_all_close


  !> Test all_zero for real and complex arrays
  subroutine test_all_zero(test_report)

    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    !> tolerance
    real(dp), parameter :: tol = 1.e-10_dp
    !> test input
    real(dp) :: a(3,4)
    !> test input
    complex(dp) :: b(3,4)

    a =  0._dp
    b = zzero

    call test_report%assert(all_zero(a,tol), 'Expected: Real array is zero')
    call test_report%assert(.not. all_zero(a + 2*tol, tol), &
                           'Expected: Real array is not zero')

    call test_report%assert(all_zero(b,tol),'Expected: Complex array is zero')
    call test_report%assert(.not. all_zero(b + 2*cmplx(2*tol,2*tol,dp), tol), &
                           'Expected: Complex array is not zero')

  end subroutine test_all_zero


  !> Test is_square
  subroutine test_is_square(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(is_square(real_symmetric_matrix_5x5), &
                           'Test is_square for real square matrix. &
                           Expected: .true.')

    call test_report%assert(.not. is_square(real_matrix_5x7), &
                           'Test is_square for real non-square matrix. &
                           Expected: .false.')

    call test_report%assert(is_square(complex_hermitian_matrix_5x5), &
                           'Test is_square for complex square matrix. &
                           Expected: .true.')

    call test_report%assert(.not. is_square(complex_matrix_5x7), &
                           'Test is_square for complex non-square matrix. &
                           Expected: .false.')
  end subroutine test_is_square


  !> Test is_hermitian
  subroutine test_is_hermitian(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(is_hermitian(real_symmetric_matrix_5x5), &
                           'Test is_square for real symmetric matrix. &
                           Expected: .true.')

    call test_report%assert(.not. is_hermitian(real_matrix_5x7(:5, :5)), &
                           'Test is_square for real square non-symmetric matrix. &
                           Expected: .false.')

    call test_report%assert(.not. is_hermitian(real_matrix_5x7), &
                           'Test is_square for real non-square matrix. &
                           Expected: .false.')

    call test_report%assert(is_hermitian(complex_hermitian_matrix_5x5), &
                           'Test is_square for complex hermitian matrix. &
                           Expected: .true.')

    call test_report%assert(.not. is_hermitian(complex_matrix_5x7(:5, :5)), &
                           'Test is_square for complex square non-hermitian matrix. &
                           Expected: .false.')

    call test_report%assert(.not. is_hermitian(complex_matrix_5x7), &
                           'Test is_square for complex non-square matrix. &
                           Expected: .false.')
  end subroutine test_is_hermitian


  !> Test is_unitary
  subroutine test_is_unitary(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(.not. is_unitary(real_symmetric_matrix_5x5, tol=10e-8_dp), &
                           'Test is_unitary for non orthogonal matrix. &
                           Expected: False')

    call test_report%assert(is_unitary(real_orthogonal_matrix_5x5, tol=10e-8_dp), &
                           'Test is_unitary for square orthogonal matrix. &
                           Expected: True')

    call test_report%assert(is_unitary(real_orthogonal_matrix_5x5(:, :3), tol=10e-8_dp), &
                           'Test is_unitary for matrix with orthogonal culumn vectors. &
                           Expected: True')

    call test_report%assert(is_unitary(real_orthogonal_matrix_5x5(2:, :), tol=10e-8_dp), &
                           'Test is_unitary for matrix with orthogonal row vectors. &
                           Expected: True')


    call test_report%assert(.not. is_unitary(complex_hermitian_matrix_5x5, tol=10e-8_dp), &
                           'Test is_unitary for non unitary matrix. &
                           Expected: False')

    call test_report%assert(is_unitary(complex_unitary_matrix_5x5, tol=10e-8_dp), &
                           'Test is_unitary for square unitary matrix. &
                           Expected: True')

    call test_report%assert(is_unitary(complex_unitary_matrix_5x5(:, :3), tol=10e-8_dp), &
                           'Test is_unitary for matrix with unitary culumn vectors. &
                           Expected: True')

    call test_report%assert(is_unitary(complex_unitary_matrix_5x5(2:, :), tol=10e-8_dp), &
                           'Test is_unitary for matrix with unitary row vectors. &
                           Expected: True')
  end subroutine test_is_unitary

  !> Test diag for complex matrix
  subroutine test_diag(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report
      
    !> tolerance
    real(dp), parameter :: tol = 1.e-10_dp
    !> Real test input
    real(dp) :: a(3,3)
    !> Complex test input
    complex(dp) :: b(3,3)

    a = 1._dp
    call test_report%assert(all_close(diag(a), [1.0_dp, 1.0_dp, 1.0_dp]), &
                      'Test diag of real matrix. &
                      & Expected: [1._dp,1._dp,1._dp] ')

    b = zone
    call test_report%assert(all_close(diag(b), [zone, zone, zone]), &
                    'Test diag of complex matrix. &
                    & Expected: [zone,zone,zone] ' )
  end subroutine test_diag


  !> Test kronecker_product for real input
  subroutine test_kronecker_product(test_report)

    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report
    !> test input
    real(dp) :: triv(2) = [1.0_dp, 1.0_dp], &
                non_triv(2) = [1.0_dp, 2.0_dp], &

                v_1(2) = [ 1.0_dp, 2.0_dp ], &
                v_2(3) = [ 3.0_dp, 4.0_dp, 5.0_dp ], &
                v_3(4) = [ 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp ], &
    !> reference
                ref(24) = [ 18.0_dp,  21.0_dp,  24.0_dp,  27.0_dp,  24.0_dp,  28.0_dp,  &
                            32.0_dp,  36.0_dp,  30.0_dp,  35.0_dp,  40.0_dp,  45.0_dp,  &
                            36.0_dp,  42.0_dp,  48.0_dp,  54.0_dp,  48.0_dp,  56.0_dp,  &
                            64.0_dp,  72.0_dp,  60.0_dp,  70.0_dp,  80.0_dp,  90.0_dp ]

    real(dp), allocatable :: kron_prod_real(:)
    complex(dp) :: ref_complex(24)
    complex(dp), allocatable :: kron_prod_complex(:)

    ! test if the function reproduces the correct structure:
    kron_prod_real = kronecker_product(triv, triv, non_triv)
    call test_report%assert(size(kron_prod_real)==8, &
                    'Test if kronecker_product is returning an array with the correct size for input vectors with the same size. &
                    Expected: 8 element vector.')

    call test_report%assert(all_close(kron_prod_real, &
                    [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp]), &
                    'Test kronecker_product structure: Third is non trivial. &
                    Expected: [1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0]')
    deallocate(kron_prod_real)

    kron_prod_real = kronecker_product(triv, non_triv, triv)
    call test_report%assert(all_close(kron_prod_real, &
                    [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]), &
                    'Test kronecker_product structure: Second input is non trivial. &
                    Expected: [1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0]')
    deallocate(kron_prod_real)

    kron_prod_real = kronecker_product(non_triv, triv, triv)
    call test_report%assert(all_close(kron_prod_real, &
                    [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]), &
                    'Test kronecker_product structure: First input is non trivial. &
                    Expected: [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]')
    deallocate(kron_prod_real)

    ! test kronecker product for real input
    kron_prod_real = kronecker_product(v_1, v_2, v_3)
    call test_report%assert(size(kron_prod_real) == 24, &
                    'Test if kronecker_product returns an array with the correct size for real input vectors with different size. &
                    Expected: 24 element vector.')

    call test_report%assert(all_close( kronecker_product(v_1, v_2, v_3), ref), &
                    'Test kronecker_product for real input: The first input has dimension 2, the second 3 and the third 4. &
                    Expected: Real vector that is the kronecker product of the input vectors.')
    deallocate(kron_prod_real)
    ! test kronecker product for complex input

    kron_prod_complex = kronecker_product(v_1*zone, v_2*zone, v_3*zone)
    ref_complex = ref * zone
    call test_report%assert(all_close(kron_prod_complex, ref_complex), &
                    'Test kronecker_product for complex input: The first input has dimension 2, the second 3 and the third 4. &
                    Expected: Complex vector that is the kronecker product of the input vectors.')
    deallocate(kron_prod_complex)

  end subroutine test_kronecker_product


  !> Test determinant
  subroutine test_determinant(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report
    !> test input
    real(dp) :: test_matrix(3, 3) = reshape([ 1.0_dp, 0.0_dp, 6.0_dp, &
                                              0.0_dp, 3.0_dp, 0.0_dp, &
                                              5.0_dp, 0.0_dp, 4.0_dp ], [ 3, 3 ])
    integer :: test_matrix_int(3, 3) 

    test_matrix_int = int(test_matrix)
           
    call test_report%assert(all_close(determinant(test_matrix), -78.0_dp), &
               'Test determinant for real bijective matrix (the deteminant is not zero). &
               Expected: -78.0')

    call test_report%assert(determinant(test_matrix_int) == -78_dp, &
               'Test determinant for integer bijective matrix (the deteminant is not zero). &
               Expected: -78')

    test_matrix = reshape([ 1.0, 2.0, 0.0, &
                            3.0, 4.0, 0.0, &
                            0.0, 0.0, 0.0 ], [ 3, 3 ])
    call test_report%assert(all_close(determinant(test_matrix), 0.0_dp), &
                            'Test determinant for real non-bijective matrix (the determinant is zero). &
                            Expected: 0.0')

    test_matrix_int =  int(test_matrix)
    call test_report%assert(determinant(test_matrix) == 0.0_dp, &
                            'Test determinant for integer non-bijective matrix (the determinant is zero). &
                            Expected: 0')

  end subroutine test_determinant


  !> Test mod1
  subroutine test_mod1(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(mod1(4, 8) ==  4, &
               'Test mod1(4, 8). &
               Expected: 4')

    call test_report%assert(mod1(5, 8) ==  5, &
               'Test mod1(5, 8). &
               Expected: 5')

    call test_report%assert(mod1(8, 4) ==  4, &
               'Test mod1(8, 4). &
               Expected: 4')

    call test_report%assert(mod1(9, 4) ==  1, &
               'Test mod1(9, 4). &
               Expected: 4')

    call test_report%assert(mod1(-4, 8) == 4, &
               'Test mod1(-4, 8). &
               Expected: 4')

    call test_report%assert(mod1(-9, 4) == -1, &
               'Test mod1(-9, 4). &
               Expected: -1')

    call test_report%assert(mod1(-7, 4) == -3, &
               'Test mod1(-7, 4). &
               Expected: -3')

    call test_report%assert(mod1(-8, 4) == 4, &
               'Test mod1(-8, 4). &
               Expected: 4')

    call test_report%assert(mod1(8, -4) == -4, &
               'Test mod1(8, -4). &
               Expected: -4')

    call test_report%assert(mod1(-5, -8) == -5, &
               'Test mod1(-5, -8). &
               Expected: -5')

    call test_report%assert(mod1(0, 5) == 5, &
               'Test mod1(0, 1). &
               Expected: 0')

    call test_report%assert(mod1(0,-3) == -3, &
               'Test mod1(0,-1). &
               Expected: 0')
  
  end subroutine test_mod1


  !> Test mask for a four element array.
  subroutine test_mask_vector(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(all_close(mask_vector([ 1.0_dp, 4.5_dp, 3.45623_dp, -5.0_dp ]*zone, [ .true., .false., .false., .true. ]), [ 1.0_dp, -5.0_dp ]*zone), &
                    'Test mask_vector for a four element array. &
                    Expected: complex two element array. The elements are the ones, corresponding to true elements in the mask ([1.0, -5.0]).')    
  end subroutine test_mask_vector

  !> Test mask for a four element array.
  subroutine test_mask(test_report)
      !> Unit test report
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(all_close(mask_vector([ 1.0_dp, 4.5_dp, 3.45623_dp, -5.0_dp ]*zone, [ .true., .false., .false., .true. ]), [ 1.0_dp, -5.0_dp ]*zone), &
                      'Test mask_vector for a four element array. &
  &                    Expected: complex two element array. The elements are the ones, corresponding to true elements in the mask: &
  &                    [1.0, -5.0]')
  end subroutine test_mask

    subroutine test_round_down(test_report)
      !> Our test object
      type(unit_test_type), intent(inout) :: test_report

      !scalar positive x; scalar positive n
      call test_report%assert(all_close(round_down(2123.77963_dp, 3), 2123.779_dp), &
                              'Test function round_down for scalar with scalar n = 3. &
                              Expected result: 2123.779')
   
      !scalar negative x; scalar positive n
      call test_report%assert(all_close(round_down(-2123.77963_dp, 3), -2123.779_dp), &
                              'Test function round_down for negative scalar with positive &
                              scalar n = 3. Expected result: -2123.779')
   
      !scalar positive x; scalar negative n
      call test_report%assert(all_close(round_down(2123.77963_dp, -2), 2100.0_dp), &
                              'Test function round_down for positive scalar with negative &
                              scalar n = -2. Expected result: 2100')
   
      !scalar negative x; scalar negative n
      call test_report%assert(all_close(round_down(-2123.77963_dp, -2), -2100.0_dp), &
                              'Test function round_down for negative scalar with negative &
                              scalar n = -2. Expected result: -2100')
      !scalar positive x; n=0
      call test_report%assert(all_close(round_down(2123.77963_dp, 0), 2123.0_dp), &
                              'Test function round_down for positive scalar with &
                              scalar n = 0. Expected result: 2123.0')
   
      !scalar negative x; n=0
      call test_report%assert(all_close(round_down(-2123.77963_dp, 0), -2123.0_dp), &
                              'Test function round_down for negative scalar with &
                              scalar n = 0. Expected result: -2123.0')
   
      !mixed array x; scalar positive n
      call test_report%assert(all_close(round_down([2123.77963_dp, -212.377963_dp, 2.12377963_dp], 3), &
                              [2123.779_dp, -212.377_dp, 2.123_dp]), 'Test function round_down for three &
                              element array with scalar n = 3. Expected result: [2123.779, -212.377, 2.123]')
   
      !positive scalar x; array positive n
      call test_report%assert(all_close(round_down(2123.77963_dp, [1, 2, 4]), [2123.7_dp, 2123.77_dp, 2123.7796_dp]), &
                              'Test function round_down for scalar with array n = [1, 2, 4]. &
                              Expected result: [2123.7, 2123.77, 2123.7796]')
   
      !array x; array n
      call test_report%assert(all_close(round_down([2123.77963_dp, 212.377963_dp, 2.12377963_dp], &
                              [2, 3, 5]), [2123.77_dp, 212.377_dp, 2.12377_dp]), &
                              'Test function round_down for array with array n = [2, 3, 5]. &
                              Expected result: [2123.77_dp, 212.377_dp, 2.12377_dp]')
   
   end subroutine test_round_down

   subroutine test_distance_matrix(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
    
    !> Test input
    !> Position A
    real(dp) :: pos_A(3, 3) = reshape([0.0_dp, 1.0_dp, 0.0_dp, &
                                      2.0_dp, 1.0_dp, 1.0_dp, &
                                      0.0_dp, 1.0_dp, 0.0_dp], [3, 3])
    !> Position B 
    real(dp) :: pos_B(3, 2) = reshape([0.0_dp, -1.0_dp, 0.0_dp,&
                                      2.0_dp, 2.0_dp, 1.0_dp], [3, 2])
    !> Position C
    real(dp) :: pos_C(3,3) = 0.0_dp    
    !> Result D
    real(dp) :: D_1(3, 2), D_2(2, 3),  D_3(3, 3)
    
    call distance_matrix(pos_A, pos_B, D_1)
    call test_report%assert(all_close(D_1, &
                            reshape([2.0_dp, 3.0_dp, 2.0_dp, &
                                    sqrt(6.0_dp), 1.0_dp, sqrt(6.0_dp)], &
                                    [3, 2])), &
                            'Tests distance_matrix with second dimension of first input matrix &
                            larger than second dimension of second input matrix dim(A) > dim(B).')

    call distance_matrix(pos_B, pos_A, D_2)
    call test_report%assert(all_close(D_2, &
                            reshape([2.0_dp, sqrt(6.0_dp), 3.0_dp, &
                                    1.0_dp, 2.0_dp, sqrt(6.0_dp)], &
                                    [2, 3])), &
                            'Tests distance_matrix with second dimension of first input matrix &
                            smaller than second dimension of second input matrix: dim(A) < dim(B).')
    
    call distance_matrix(pos_A, pos_C, D_3)    
    call test_report%assert(all_close(D_3, &
                            reshape([1.0_dp, sqrt(6.0_dp), 1.0_dp, &
                                    1.0_dp, sqrt(6.0_dp), 1.0_dp, &
                                    1.0_dp, sqrt(6.0_dp), 1.0_dp], &
                                    [3, 3])), &
                            'Tests distance_matrix with second dimension of first input matrix &
                            equal to second dimension of second input matrix: dim(A) = dim(B).')

  end subroutine test_distance_matrix

  subroutine test_outer_sum(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: C_1(3, 3), C_2(3, 2), C_3(2, 3)

    call outer_sum([-1.0_dp, 2.0_dp, 5.0_dp], [1.0_dp, 2.0_dp, 3.0_dp], C_1)
    call test_report%assert(all_close(C_1, reshape([0.0_dp, 3.0_dp, 6.0_dp, 1.0_dp, &
                            4.0_dp, 7.0_dp, 2.0_dp, 5.0_dp, 8.0_dp], [3,3])), &
                            'Test function outer_sum for array size a == b.')

    call outer_sum([1.0_dp, 2.0_dp, 3.0_dp],[-1.0_dp, 2.0_dp], C_2)
    call test_report%assert(all_close(C_2, reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [3,2])), &
                    'Test function outer_sum for array size a > b.')
    
    call outer_sum([-1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp, 3.0_dp], C_3)
    call test_report%assert(all_close(C_3, reshape([0.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 5.0_dp], [2,3])), &
                           'Test function outer_sum for array size a < b.')

  end subroutine test_outer_sum

  subroutine test_get_subinterval_indices(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: res1(2, 2), res2(2, 1)
    integer :: vector_length, block_size
    
    vector_length = 6
    block_size = 3
    res1 = get_subinterval_indices(vector_length, block_size)
    call test_report%assert(size(res1, 1) == 2 .and. size(res1, 2) == 2, &
                            'Tests returned size of function get_subinterval_indices.')
    call test_report%assert(all_close(res1(:, 1), [1.0_dp, 3.0_dp]), "First set of limits.")
    call test_report%assert(all_close(res1(:, 2), [4.0_dp, 6.0_dp]), "Second set of limits.")
    
    vector_length = 6
    block_size = 6
    res2 = get_subinterval_indices(vector_length, block_size)
    call test_report%assert(size(res2, 1) == 2 .and. size(res2, 2) == 1, &
                            'Tests returned size of function get_subinterval_indices for an equal block size.')
    call test_report%assert(all_close(res2(:, 1), [1.0_dp, 6.0_dp]), &
                            "Set of limits for function get_subinterval_indices for an equal block size.")
                 
  end subroutine test_get_subinterval_indices

end module math_utils_test
