!> Unit tests for [[seed_generation]].
module seed_generation_test
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use precision, only: dp
  use seed_generation, only: set_seed

  implicit none

  private
  public :: seed_generation_test_driver

  contains

  !> Run tests for seed generation
  subroutine seed_generation_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 4

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_seed_int(test_report)
    call test_seed_wrapper(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report("linear_algebra_3d", kill_on_failure)
    else
      call test_report%report("linear_algebra_3d")
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine seed_generation_test_driver

  subroutine test_seed_int(test_report)
    !> Test report
    type(unit_test_type), intent(inout)  :: test_report

    real(dp) :: random_number_1, random_number_2

    call set_seed(1345)
    call random_number(random_number_1)

    call set_seed(1345)
    call random_number(random_number_2)

    call test_report%assert(random_number_1 == random_number_2, &
            "random_number fixed when using the same seeds.")

    call set_seed(1345)
    call random_number(random_number_1)

    call set_seed(1343525)
    call random_number(random_number_2)

    call test_report%assert(random_number_1 /= random_number_2, &
            "random_number not fixed when using different seeds.")
  end subroutine

  subroutine test_seed_wrapper(test_report)
    !> Test report
    type(unit_test_type), intent(inout)  :: test_report

    real(dp) :: random_number_1, random_number_2

    call set_seed("fixed")
    call random_number(random_number_1)

    call set_seed("fixed")
    call random_number(random_number_2)

    call test_report%assert(random_number_1 == random_number_2, &
            "random_number fixed when using fixed seed.")

    call set_seed("clock")
    call random_number(random_number_1)
    call sleep(1)

    call set_seed("clock")
    call random_number(random_number_2)

    call test_report%assert(random_number_1 /= random_number_2, &
            "random_number not fixed when using system time as seed.")
  end subroutine test_seed_wrapper

end module seed_generation_test