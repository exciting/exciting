!> Module for collecting unit test drivers for modules in the math directory.
module math_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here
  use math_utils_test, only: math_utils_test_driver
  use grid_utils_test, only: grid_utils_test_driver
  use linear_algebra_3d_test, only: linear_algebra_3d_test_driver
  use matrix_exp_test, only: matrix_exp_test_driver
  use integration_test, only: integration_test_driver
  use normalize_test, only: normalize_test_driver
  use sorting_test, only: sorting_test_driver
  use multi_index_conversion_test, only: multi_index_conversion_test_driver
  use seed_generation_test, only: seed_generation_test_driver
  use sh_product_test, only: sh_product_test_driver
  use projection_test, only: projection_test_driver
  use matrix_contraction_test, only: matrix_contraction_test_driver
  use iterative_solver_test, only: iterative_solver_test_driver

  private
  public :: math_test_driver

  contains

  subroutine math_test_driver(mpiglobal, kill_on_failure_)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure_ 

    logical :: kill_on_failure
    logical, parameter :: kill_on_failure_default = .true.

    kill_on_failure = kill_on_failure_default
    if( present(kill_on_failure_) ) kill_on_failure = kill_on_failure_

    ! Call test drivers here
    call math_utils_test_driver(mpiglobal, kill_on_failure)
    call grid_utils_test_driver(mpiglobal, kill_on_failure)
    call linear_algebra_3d_test_driver(mpiglobal, kill_on_failure)
    call matrix_exp_test_driver(mpiglobal, kill_on_failure)
    call integration_test_driver(mpiglobal, kill_on_failure)
    call normalize_test_driver(mpiglobal, kill_on_failure)
    call sorting_test_driver(mpiglobal, kill_on_failure)
    call multi_index_conversion_test_driver(mpiglobal, kill_on_failure)
    call seed_generation_test_driver(mpiglobal, kill_on_failure)
    call sh_product_test_driver(mpiglobal, kill_on_failure)
    call projection_test_driver(mpiglobal, kill_on_failure)
    call matrix_contraction_test_driver(mpiglobal, kill_on_failure)
    call iterative_solver_test_driver(mpiglobal, kill_on_failure)
  end subroutine math_test_driver

end module math_test_drivers
