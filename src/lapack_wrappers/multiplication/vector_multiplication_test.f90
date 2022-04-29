module vector_multiplication_test
  use precision, only: dp
  use constants, only: zi
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close
  use mock_arrays, only: real_vector_5, &
                         complex_vector_5, &
                         real_vector_7, &
                         complex_vector_7, &
                         real_matrix_5x7, &
                         complex_matrix_5x7
  use vector_multiplication, only: norm, dot_multiply, outer_product

  implicit none

  private
  public :: vector_multiplication_test_driver
  
contains

  !> Run tests for vector multiplication
  subroutine vector_multiplication_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 14

    call test_report%init(n_assertions, mpiglobal)

    ! Run unit tests

    call test_dot_multiply(test_report)

    call test_outer_product(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('vector_multiplication', kill_on_failure)
    else
      call test_report%report('vector_multiplication')
    end if

    call test_report%finalise()
  end subroutine vector_multiplication_test_driver


  !> Test norm
  subroutine test_norm(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report

    real(dp) :: norm_ref
    
    real(dp), allocatable :: p(:)
    complex(dp), allocatable :: phi(:)

    p = real_vector_5
    phi = complex_vector_5

    norm_ref = sqrt(sum(p ** 2))
    call test_report%assert(all_close(norm(p), norm_ref), &
                           'Test norm for a real vector.')
    
    norm_ref = sqrt(sum(real(conjg(phi) * phi, dp)))                           
    call test_report%assert(all_close(norm(phi), norm_ref), &
                           'Test norm for a complex vector.')
  end subroutine test_norm


  !> Test dot product
  subroutine test_dot_multiply(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report

    real(dp), allocatable :: p(:), q(:)
    complex(dp), allocatable :: phi(:), psi(:)

    p = real_vector_5
    q = real_vector_7(1 : 5)
    phi = complex_vector_5
    psi = complex_vector_7(1 : 5)

    call test_report%assert(all_close(dot_multiply(p, q), dot_multiply(p, q)), &
                           'Test dot_multiply real vectors. &
                           Expeceted: p * q')
    
    call test_report%assert(all_close(dot_multiply(phi, psi), dot_multiply(phi, psi)), &
                           'Test dot_multiply for complex vectors. &
                           Expeceted: conjg(phi) * psi')

    ! dot_multiply conjugates the left vector by default. 
    call test_report%assert(all_close(dot_multiply(phi, psi, conjg_a = .false.), dot_multiply(conjg(phi), psi)), & 
                           'Test dot_multiply for complex vectors. &
                           Expeceted: phi * psi')

    call test_report%assert(all_close(dot_multiply(p, psi), dot_multiply(p, psi)), &
                           'Test dot_multiply for real times complex vectors. &
                           Expeceted: phi * psi')

    call test_report%assert(all_close(dot_multiply(phi, q), dot_multiply(phi, q)), &
                           'Test dot_multiply for complex times real vectors. &
                           Expeceted: phi * psi')

    ! dot_multiply conjugates the left side by default. 
    call test_report%assert(all_close(dot_multiply(phi, q, conjg_a = .false.), dot_multiply(conjg(phi), q)), &
                           'Test dot_multiply for complex times real vectors. &
                           Expeceted: phi * psi')
  end subroutine test_dot_multiply


  !> Test outer product
  subroutine test_outer_product(test_report)
    !> Test object
    type(unit_test_type), intent(inout) :: test_report

    integer :: i, j
    real(dp) :: p(5), q(7), result_real(5, 7), reference_real(5, 7)
    complex(dp) :: phi(5), psi(7), result_complex(5, 7), reference_complex(5, 7)

    p = real_vector_5
    q = real_vector_7
    result_real = real_matrix_5x7
    reference_real = real_matrix_5x7

    phi = complex_vector_5
    psi = complex_vector_7
    result_complex = complex_matrix_5x7
    reference_complex = complex_matrix_5x7

    ! real-real
    do i=1, size(p)
      do j=1, size(q)
        reference_real(i, j) = reference_real(i, j) + p(i) * q(j)
      end do 
    end do
    call outer_product(p, q, result_real)
    call test_report%assert(all_close(result_real, reference_real), &
                           'Test outer product subroutine call for real vectors. &
                           Expected: C + p * transpose(q)')
          
    ! complex-complex
    do i=1, size(phi)
      do j=1, size(psi)
        reference_complex(i, j) = reference_complex(i, j) + phi(i) * psi(j)
      end do 
    end do
    call  outer_product(phi, psi, result_complex)
    call test_report%assert(all_close(result_complex, reference_complex), &
                           'Test outer product subroutine call for complex vectors. &
                           Expeceted: C + phi * transpose(psi)')   
    
    ! complex-complex conjugate
    do i=1, size(phi)
      do j=1, size(psi)
        reference_complex(i, j) = reference_complex(i, j) + phi(i) * conjg(psi(j))
      end do 
    end do
    call  outer_product(phi, psi, result_complex, conjg_b = .true.)
    call test_report%assert(all_close(result_complex, reference_complex), &
                           'Test outer product subroutine call for complex vectors. &
                           Expeceted: C + phi * transpose(psi)')

    ! real-complex
    do i=1, size(p)
      do j=1, size(psi)
        reference_complex(i, j) = reference_complex(i, j) + p(i) * psi(j)
      end do 
    end do
    call outer_product(p, psi, result_complex)
    call test_report%assert(all_close(result_complex, reference_complex), &
                           'Test outer product subroutine call for a real times a complex vector. &
                           Expeceted: p * transpose(psi)')

    ! real-complex conjugate
    do i=1, size(p)
      do j=1, size(psi)
        reference_complex(i, j) = reference_complex(i, j) + p(i) * conjg(psi(j))
      end do 
    end do
    call outer_product(p, psi, result_complex, conjg_b = .true.)
    call test_report%assert(all_close(result_complex, reference_complex), &
                           'Test outer product subroutine call for a real times a complex vector. &
                           Expeceted: p * adjungate(psi)')                       

    ! complex-real
    do i=1, size(phi)
      do j=1, size(q)
        reference_complex(i, j) = reference_complex(i, j) + phi(i) * q(j)
      end do 
    end do
    call outer_product(phi, q, result_complex)
    call test_report%assert(all_close(result_complex, reference_complex), &
                           'Test outer product for a complex times a real vector. &
                           Expeceted: phi * transpose(q)')
  end subroutine test_outer_product

end module vector_multiplication_test
