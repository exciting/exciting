!> Unit tests for singular_value_decomposition
module singular_value_decomposition_test
  use precision, only: dp
  use constants, only: zone, zzero
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, is_unitary    
  use mock_arrays, only: real_matrix_5x7, &
                         real_matrix_7x5, &
                         real_full_rank_matrix_5x5, &
                         real_rank_4_matrix_5x7, &
                         real_rank_4_matrix_7x5, &
                         complex_matrix_5x7, &
                         complex_matrix_7x5, &
                         complex_full_rank_matrix_5x5, &
                         complex_rank_4_matrix_5x7, &
                         complex_rank_4_matrix_7x5
  use singular_value_decomposition, only: svd_divide_conquer

  implicit none

  private
  public :: singular_value_decomposition_test_driver
  
contains

  !> Run tests for singular value decomposition.
  subroutine singular_value_decomposition_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 36

    call test_report%init(n_assertions, mpiglobal)

    call test_singular_value_decomposition_real(test_report)
    
    call test_singular_value_decomposition_complex(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('singular_value_decomposition', kill_on_failure)
    else
      call test_report%report('singular_value_decomposition')
    end if

    call test_report%finalise()
  end subroutine singular_value_decomposition_test_driver


  !> Test singular value decomposition for real input.
  subroutine test_singular_value_decomposition_real(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    integer :: i
    real(dp), allocatable :: sigma(:), sigma_ref(:), A(:, :), U(:, :), V_T(:, :), &
                             sigma_matrix(:, :), U_sigma_V_T(:, :), &
                             U_ref(:, :), V_T_ref(:, :)

    A = real_matrix_5x7
    
    ! Calculate all left and right singular vectors
    allocate(sigma(5))
    allocate(U(5, 5))
    allocate(V_T(7, 7))
    call svd_divide_conquer(A, sigma, U, V_T)
 
    U_ref = U
    V_T_ref = V_T
    
    ! U * U_T = 1
    call test_report%assert(is_unitary(U), &
                              'Test for svd_divide_conquer with real input and calculating all singular vectors &
                              that the left singular vectors are an orthorgonal matrix. &
                              Expected: U * U**T = I')

    ! V * V_T = 1
    call test_report%assert(is_unitary(V_T), &
                              'Test for svd_divide_conquer with real input and calculating all singular vectors &
                              that the right singular vectors are an orthogonal matrix. &
                              Expected: V * V**T = I')

    ! U * sigma * V_T = A
    allocate(sigma_matrix(5, 7))
    sigma_matrix = zzero
    do i = 1, 5
      sigma_matrix(i, i) = zone * sigma(i)
    end do
    
    U_sigma_V_T = matmul(U, matmul(sigma_matrix, V_T))
    call test_report%assert(all_close(A, U_sigma_V_T), &
                           'Test for svd_divide_conquer with complex input and calculating all singular vectors &
                           that the definition is fullfilled. &
                           Expected: A = U * sigma * V**T')


    ! Calculate only the first 5 right singular vectors
    deallocate(V_T); allocate(V_T(5, 7))
    call svd_divide_conquer(A, sigma, U, V_T)

    ! U * U_T = 1
    call test_report%assert(is_unitary(U), &
                              'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                              right singular vectors that the left singular vectors are an orthorgonal matrix. &
                              Expected: U * U**T = I')

    ! V * V_T = 1
    call test_report%assert(is_unitary(V_T), &
                              'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                              right singular vectors that the right singular vectors are an orthogonal matrix. &
                              Expected: V**T * V = I')

    ! U = U_ref
    call test_report%assert(all_close(U, U_ref), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           right singular vectors that the left singular vectors are correct. &
                           Expected:  U = U_ref')

    ! V_T = V_T_ref(:, :5)
    call test_report%assert(all_close(V_T, V_T_ref(:5, :)), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           right singular vectors that the right singular vectors are correct. &
                           Expected: V_T = V_T_ref(:5, :)')

    ! U * sigma * V_T = A
    deallocate(sigma_matrix); allocate(sigma_matrix(5, 5))
    sigma_matrix = zzero
    do i=1, 5
      sigma_matrix(i, i) = zone * sigma(i)
    end do
    
    U_sigma_V_T = matmul(U, matmul(sigma_matrix, V_T))
    call test_report%assert(all_close(A, U_sigma_V_T), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           right singular vectors that the definition is fullfilled. &
                           Expected: A = U * sigma * V**T')


    ! Calculate only the first 5 left singular vectors
    deallocate(U, V_T); allocate(U(7, 7)); allocate(V_T(5, 5)) 
    A = real_matrix_7x5
    call svd_divide_conquer(A, sigma, U, V_T)
    U_ref = U; V_T_ref = V_T 

    deallocate(U); allocate(U(7, 5))
    call svd_divide_conquer(A, sigma, U, V_T)

    ! U * U_T = 1
    call test_report%assert(is_unitary(U), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           left singular vectors that the left singular vectors are an orthorgonal matrix. &
                           Expected: U**T * U = I')

    ! V * V_T = 1
    call test_report%assert(is_unitary(V_T), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           left singular vectors that the right singular vectors are an orthogonal matrix. &
                           Expected: V * V**T = I')

    ! U = U_ref
    call test_report%assert(all_close(U, U_ref(:, :5)), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           left singular vectors that the left singular vectors are correct. &
                           Expected:  U = U_ref(:, :)')

    ! V_T = V_T_ref(:5, :)
    call test_report%assert(all_close(V_T, V_T_ref), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           left singular vectors that the right singular vectors are correct. &
                           Expected: V_T = V_T_ref(:5, :)')

    ! U * sigma * V_T = A
    deallocate(sigma_matrix); allocate(sigma_matrix(5, 5))
    sigma_matrix = zzero
    do i=1, 5
      sigma_matrix(i, i) = zone * sigma(i)
    end do
    
    U_sigma_V_T = matmul(U, matmul(sigma_matrix, V_T))
    call test_report%assert(all_close(A, U_sigma_V_T), &
                           'Test for svd_divide_conquer with real input and calculating only the first 5 &
                           left singular vectors that the definition is fullfilled. &
                           Expected: A = U * sigma * V**T')

    ! Calculate only sigma
    sigma_ref = sigma 
    call svd_divide_conquer(A, sigma)

    call test_report%assert(all_close(sigma, sigma_ref), &
                           'Test svd_divide_conquer with real input and not calculating any singular vectors. &
                           Expected: Sigma should be the same as before.')


  end subroutine test_singular_value_decomposition_real

  
  !> Test singular value decomposition for complex input.
  subroutine test_singular_value_decomposition_complex(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    integer :: i
    real(dp), allocatable :: sigma(:), sigma_ref(:)
    complex(dp), allocatable :: A(:, :), U(:, :), V_dagger(:, :), &
                                sigma_matrix(:, :), U_sigma_V_dagger(:, :), &
                                U_ref(:, :), V_dagger_ref(:, :)

    A = complex_matrix_5x7
    
    ! Calculate all left and right singular vectors
    allocate(sigma(5))
    allocate(U(5, 5))
    allocate(V_dagger(7, 7))
    call svd_divide_conquer(A, sigma, U, V_dagger)

    U_ref = U 
    V_dagger_ref = V_dagger
    
    ! U * U_dagger = 1
    call test_report%assert(is_unitary(U), &
                              'Test for svd_divide_conquer with complex input that the left singular vectors &
                              are an unitary matrix. &
                              Expected: U * U**C = I')

    ! V * V_dagger = 1             
    call test_report%assert(is_unitary(V_dagger), &
                              'Test for complex svd_divide_conquer that the right singular vectors &
                              are an unitary matrix. &
                              Expected: V * V**C = I')

    ! U * sigma * V_dagger = A
    allocate(sigma_matrix(5, 7))
    sigma_matrix = zzero
    do i=1, 5
      sigma_matrix(i, i) = zone * sigma(i)
    end do
    
    U_sigma_V_dagger = matmul(U, matmul(sigma_matrix, V_dagger))
    call test_report%assert(all_close(A, U_sigma_V_dagger), &
                           'Test for svd_divide_conquerv with complex input that the definition is fullfilled. &
                           Expected: A = U * sigma * V**C')

    
    ! Calculate only the first 5 right singular vectors
    deallocate(V_dagger); allocate(V_dagger(5, 7))
    call svd_divide_conquer(A, sigma, U, V_dagger)

    ! U * U_dagger = 1
    call test_report%assert(is_unitary(U), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           right singular vectors that the left singular vectors are an unitary matrix. &
                           Expected: U * U**C = I')

    ! V * V_dagger = 1
    call test_report%assert(is_unitary(V_dagger), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           right singular vectors that the right singular vectors are an unitary matrix. &
                           Expected: V**C * V = I')

    ! U = U_ref
    call test_report%assert(all_close(U, U_ref), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           right singular vectors that the right singular vectors are correct. &
                           Expected:  U = U_ref')

    ! V_dagger = V_dagger_ref(:, :5)
    call test_report%assert(all_close(V_dagger, V_dagger_ref(:5, :)), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           right singular vectors that the right singular vectors are correct. &
                           Expected: V_dagger = V_dagger_ref(:, :5)')

    ! U * sigma * V_dagger = A
    deallocate(sigma_matrix); allocate(sigma_matrix(5, 5))
    sigma_matrix = zzero
    do i=1, 5
      sigma_matrix(i, i) = zone * sigma(i)
    end do
    
    U_sigma_V_dagger = matmul(U, matmul(sigma_matrix, V_dagger))
    call test_report%assert(all_close(A, U_sigma_V_dagger), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           right singular vectors that the definition is fullfilled. &
                           Expected: A = U * sigma * V**C')


    ! Calculate only the first 5 left singular vectors
    deallocate(U, V_dagger); allocate(U(7, 7)); allocate(V_dagger(5, 5)) 
    A = real_matrix_7x5
    call svd_divide_conquer(A, sigma, U, V_dagger)

    U_ref = U
    V_dagger_ref = V_dagger

    deallocate(U); allocate(U(7, 5))
    call svd_divide_conquer(A, sigma, U, V_dagger)

    ! U * U_dagger = 1
    call test_report%assert(is_unitary(U), &
                              'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                              left singular vectors that the left singular vectors are an unitary matrix. &
                              Expected: U**C * U = I')

    ! V * V_dagger = 1
    call test_report%assert(is_unitary(V_dagger), &
                              'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                              left singular vectors that the right singular vectors are an unitary matrix. &
                              Expected: V * V**C = I')

    ! U = U_ref(:, :5)
    call test_report%assert(all_close(U, U_ref(:, :5)), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           left singular vectors that the right singular vectors are correct. &
                           Expected:  U = U_ref(:, :5)')

    ! V_dagger = V_dagger_ref
    call test_report%assert(all_close(V_dagger, V_dagger_ref), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           left singular vectors that the right singular vectors are correct. &
                           Expected: V_dagger = V_dagger_ref')

    ! U * sigma * V_dagger = A
    deallocate(sigma_matrix); allocate(sigma_matrix(5, 5))
    sigma_matrix = zzero
    do i=1, 5
      sigma_matrix(i, i) = zone * sigma(i)
    end do
    
    U_sigma_V_dagger = matmul(U, matmul(sigma_matrix, V_dagger))
    call test_report%assert(all_close(A, U_sigma_V_dagger), &
                           'Test for svd_divide_conquer with complex input and calculating only the first 5 &
                           left singular vectors that the definition is fullfilled. &
                           Expected: A = U * sigma * V**C')

    ! Calculate only sigma
    sigma_ref = sigma 
    call svd_divide_conquer(A, sigma)

    call test_report%assert(all_close(sigma, sigma_ref), &
                           'Test svd_divide_conquer with complex input and not calculating any singular vectors. &
                           Expected: Sigma should be the same as before.')
  end subroutine test_singular_value_decomposition_complex

end module singular_value_decomposition_test