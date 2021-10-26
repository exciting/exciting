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
                        kronecker_product, &
                        determinant, &
                        permanent, &
                        mod1, &
                        shuffle_vector, &
                        mask_vector
  use mock_arrays, only: real_matrix_5x7, &
                         real_symmetric_matrix_5x5, &
                         complex_hermitian_matrix_5x5, &
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
    integer, parameter :: n_assertions = 56

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests

    call test_all_close(test_report)

    call test_all_zero(test_report)

    call test_is_square(test_report)

    call test_is_hermitian(test_report)

    call test_diagonal(test_report)
    
    call test_kronecker_product(test_report)
    
    call test_determinant(test_report)
    
    call test_mod1(test_report)

    call test_mask_vector(test_report)

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

    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    !> tolerance
    real(dp), parameter :: tol = 1.e-10_dp
    !> Real test input
    real(dp) :: a(3,4) 
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

    !> Our test object
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
    !> Our test object
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
    !> Our test object
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


  !> Test diagonal for complex matrix
  subroutine test_diagonal(test_report)

    !> Our test object
    type(unit_test_type), intent(inout) :: test_report
      
    !> tolerance
    real(dp), parameter :: tol = 1.e-10_dp
    !> Real test input
    real(dp) :: a(3,3)
    !> Complex test input
    complex(dp) :: b(3,3)

    a = 1._dp
    call test_report%assert(all_close(diag(a), [1.0_dp, 1.0_dp, 1.0_dp]), &
                      'Test diagonal of real matrix. &
                      & Expected: [1._dp,1._dp,1._dp] ')

    b = zone
    call test_report%assert(all_close(diag(b), [zone, zone, zone]), &
                    'Test diagonal of complex matrix. &
                    & Expected: [zone,zone,zone] ' )
  end subroutine test_diagonal

  
  !> Test kronecker_product for real input
  subroutine test_kronecker_product(test_report)

    !> Our test object
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
    !> Our test object
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

  subroutine test_mod1(test_report)
    !> Our test object
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
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(all_close(mask_vector([ 1.0_dp, 4.5_dp, 3.45623_dp, -5.0_dp ]*zone, [ .true., .false., .false., .true. ]), [ 1.0_dp, -5.0_dp ]*zone), &
                    'Test mask_vector for a four element array. &
                    Expected: complex two element array. The elements are the ones, corresponding to true elements in the mask ([1.0, -5.0]).')    
  end subroutine test_mask_vector

  !> Test mask for a four element array.
  subroutine test_mask(test_report)
      !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    call test_report%assert(all_close(mask_vector([ 1.0_dp, 4.5_dp, 3.45623_dp, -5.0_dp ]*zone, [ .true., .false., .false., .true. ]), [ 1.0_dp, -5.0_dp ]*zone), &
                      'Test mask_vector for a four element array. &
  &                    Expected: complex two element array. The elements are the ones, corresponding to true elements in the mask: &
  &                    [1.0, -5.0]')
  end subroutine test_mask

end module math_utils_test
