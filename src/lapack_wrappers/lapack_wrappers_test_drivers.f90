!> Module for collecting unit test drivers for modules in the lapack_wrappers directory.
module lapack_wrappers_test_drivers
  use modmpi, only: mpiinfo
  ! multiplication
  use vector_multiplication_test, only: vector_multiplication_test_driver
  use general_matrix_multiplication_test, only: general_matrix_multiplication_test_driver
  use hermitian_matrix_multiplication_test, only: hermitian_matrix_multiplication_test_driver
  ! decomposition
  use singular_value_decomposition_test, only: singular_value_decomposition_test_driver
  use lu_factorization_test, only: lu_factorization_test_driver
  use qr_factorization_test, only: qr_factorization_test_driver
  ! diagonalization
  use diagonalize_tridiagonal_test, only: diagonalize_tridiagonal_test_driver
  ! utils
  use matrix_rank_test, only: matrix_rank_test_driver
  use determinant_test, only: determinant_test_driver
  use inverse_test, only: inverse_test_driver
  ! Tensor contraction
  use tensor_contraction_test, only: tensor_contraction_test_driver

  implicit none
  
  private
  public :: lapack_wrappers_test_driver

contains

  subroutine lapack_wrappers_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure 

    ! multiplication
    call vector_multiplication_test_driver(mpiglobal, kill_on_failure)
    call general_matrix_multiplication_test_driver(mpiglobal, kill_on_failure)
    call hermitian_matrix_multiplication_test_driver(mpiglobal, kill_on_failure)
    ! decomposition
    call singular_value_decomposition_test_driver(mpiglobal, kill_on_failure)
    call lu_factorization_test_driver(mpiglobal, kill_on_failure)
    call qr_factorization_test_driver(mpiglobal, kill_on_failure)
    ! Diagonalization wrappers
    call diagonalize_tridiagonal_test_driver(mpiglobal, kill_on_failure)
    ! utils
    call matrix_rank_test_driver(mpiglobal, kill_on_failure)
    call determinant_test_driver(mpiglobal, kill_on_failure)
    call inverse_test_driver(mpiglobal, kill_on_failure)
    call tensor_contraction_test_driver(mpiglobal, kill_on_failure)

  end subroutine lapack_wrappers_test_driver

end module lapack_wrappers_test_drivers