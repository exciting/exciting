module vector_multiplication_test
  use precision, only: dp
  use constants, only: zi
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close
  use mock_arrays, only: real_vector_5, &
                         complex_vector_5, &
                         real_vector_7, &
                         complex_vector_7
  use vector_multiplication, only: dot_multiply

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

    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 6

    call test_report%init(n_assertions, mpiglobal)

    ! Run unit tests

    call test_dot_multiply(test_report)

    

    if (present(kill_on_failure)) then
      call test_report%report('vector_multiplication', kill_on_failure)
    else
      call test_report%report('vector_multiplication')
    end if

    call test_report%finalise()
  end subroutine vector_multiplication_test_driver


  subroutine test_dot_multiply(test_report)
    !> Test object
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

end module vector_multiplication_test