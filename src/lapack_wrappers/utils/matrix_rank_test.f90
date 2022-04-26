module matrix_rank_test
  use precision, only: dp
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use mock_arrays, only: real_full_rank_matrix_5x5, &
                         real_rank_4_matrix_5x7, &
                         real_rank_4_matrix_7x5, &
                         complex_full_rank_matrix_5x5, &
                         complex_rank_4_matrix_5x7, &
                         complex_rank_4_matrix_7x5
  use matrix_rank, only: matrix_rank_SVD

  implicit none

  private
  public :: matrix_rank_test_driver

contains

  !> Run tests for singular value decomposition.
  subroutine matrix_rank_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 8

    call test_report%init(n_assertions, mpiglobal)

    call test_matrix_rank_SVD(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('matrix_rank_lapack', kill_on_failure)
    else
      call test_report%report('matrix_rank_lapack')
    end if

    call test_report%finalise()
  end subroutine matrix_rank_test_driver


  !> Test matrix rank relying on singular value decomposition.
  subroutine test_matrix_rank_SVD(test_report)
    !> Test report object
    type(unit_test_type) :: test_report

    call test_report%assert(matrix_rank_SVD(real_full_rank_matrix_5x5) == 5, &
                           'Test matrix_rank_SVD for a real full rank matrix. &
                           Expected: The rank is 5.')

    call test_report%assert(matrix_rank_SVD(real_rank_4_matrix_5x7) == 4, &
                           'Test matrix_rank_SVD for real 5 x 7 matrix with rank 4. &
                           Expected: The rank is 4.')

    call test_report%assert(matrix_rank_SVD(real_rank_4_matrix_7x5) == 4, &
                           'Test matrix_rank_SVD for real 7 x 5 matrix with rank 4. &
                           Expected: The rank is 4.')

    call test_report%assert(matrix_rank_SVD(real_full_rank_matrix_5x5, tol=1._dp) == 4, &
                           'Test matrix_rank_SVD for a real full rank matrix with a given tolerance for &
                           defining the singular values as zero. &
                           Expected: The rank is 4.')

    call test_report%assert(matrix_rank_SVD(complex_full_rank_matrix_5x5) == 5, &
                           'Test matrix_rank_SVD for a complex full rank matrix. &
                           Expected: The rank is 5.')

    call test_report%assert(matrix_rank_SVD(complex_rank_4_matrix_5x7) == 4, &
                           'Test matrix_rank_SVD for complex 5 x 7 matrix with rank 4. &
                           Expected: The rank is 4.')

    call test_report%assert(matrix_rank_SVD(complex_rank_4_matrix_7x5) == 4, &
                           'Test matrix_rank_SVD for complex 7 x 5 matrix with rank 4. &
                           Expected: The rank is 4.')

    call test_report%assert(matrix_rank_SVD(complex_full_rank_matrix_5x5, tol=10._dp) == 4, &
                           'Test matrix_rank_SVD for a complex full rank matrix with a given tolerance for &
                           defining the singular values as zero. &
                           Expected: The rank is 4.')
  end subroutine test_matrix_rank_SVD

end module matrix_rank_test