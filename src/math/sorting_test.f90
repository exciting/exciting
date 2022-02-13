!> Unit tests for sorting.f90
module sorting_test
   use precision, only: dp
   use modmpi, only: mpiinfo
   use unit_test_framework, only : unit_test_type
   use sorting, only: sort_index_1d, sort_index_2d

   implicit none
   private

   public :: sorting_test_driver

contains

  !> Driver for sorting unit tests
  subroutine sorting_test_driver(mpiglobal, kill_on_failure)
    type(mpiinfo), intent(in) :: mpiglobal
    logical, intent(in), optional :: kill_on_failure

    type(unit_test_type) :: test_report
    logical :: kill
    integer, parameter :: n_assertions = 12

    if (present(kill_on_failure)) then
        kill = merge(.true., .false., kill_on_failure)
    end if

    call test_report%init(n_assertions, mpiglobal)

    call test_sort_index_1d_integer(test_report)
    call test_sort_index_1d_realdp(test_report) 
    call test_sort_index_2d(test_report) 
    
    call test_report%report('sorting', kill)
    call test_report%finalise()

  end subroutine 


  !> Test sorting for 1D array of integers
  subroutine test_sort_index_1d_integer(test_report)
    type(unit_test_type), intent(inout) :: test_report

    !> Unordered vector of integers
    integer :: a(6)
    !> Results
    integer :: indices(6), a_sorted(6)
    logical :: valid_sorting

    ! Unordered, positive integers
    a = [2, 5, 6, 3, 0, 1]
    indices = sort_index_1d(size(a), a)
    a_sorted = a(indices)
    call test_report%assert(all(indices == [5, 6, 1, 4, 2, 3]), message="Indices for array reordering")
    call test_report%assert(all(a_sorted == [0, 1, 2, 3, 5, 6]), &
        message="Expect array ordered lowest to highest")

    ! Unordered, negative integers
    a = [1, 2, 0, -2, -1, -3]
    indices = sort_index_1d(size(a), a)
    a_sorted = a(indices)
    call test_report%assert(all(a_sorted == [-3, -2, -1, 0, 1, 2]), &
        message="Expect array ordered lowest to highest")

    ! Unordered integers with duplicates
    a = [0, 2, 0, -2, -1, 1]
    indices = sort_index_1d(size(a), a)
    valid_sorting = all(indices == [4, 5, 1, 3, 6, 2]) .or. all(indices == [4, 5, 3, 1, 6, 2])
    call test_report%assert(valid_sorting, message="Degenerate values can be returned in any valid order")

    a_sorted = a(indices)
    call test_report%assert(all(a_sorted == [-2, -1, 0, 0, 1, 2]), &
        message="Expect array ordered lowest to highest")

  end subroutine
  
  !> Test sorting for 1D array of reals
  subroutine test_sort_index_1d_realdp(test_report)
    type(unit_test_type), intent(inout) :: test_report

    !> Unordered vector
    real(dp) :: a(6)
    !> Results
    integer :: indices(6)
    real(dp) :: a_sorted(6)
    logical :: valid_sorting

    ! Unordered, positive reals
    a = [2._dp, 5._dp, 6._dp, 3._dp, 0._dp, 1._dp]
    indices = sort_index_1d(size(a), a)
    a_sorted = a(indices)
    call test_report%assert(all(indices == [5, 6, 1, 4, 2, 3]), message="Indices for array reordering")
    call test_report%assert(all(a_sorted == [0._dp, 1._dp, 2._dp, 3._dp, 5._dp, 6._dp]), &
        message="Expect array ordered lowest to highest")

    ! Unordered, negative reals
    a = [1._dp, 2._dp, 0._dp, -2._dp, -1._dp, -3._dp]
    indices = sort_index_1d(size(a), a)
    a_sorted = a(indices)
    call test_report%assert(all(a_sorted == [-3._dp, -2._dp, -1._dp, 0._dp, 1._dp, 2._dp]), &
        message="Expect array ordered lowest to highest")

    ! Unordered reals with duplicates
    a = [0._dp, 2._dp, 0._dp, -2._dp, -1._dp, 1._dp]
    indices = sort_index_1d(size(a), a)
    valid_sorting = all(indices == [4, 5, 1, 3, 6, 2]) .or. all(indices == [4, 5, 3, 1, 6, 2])
    call test_report%assert(valid_sorting, message="Degenerate values can be returned in any valid order")
    a_sorted = a(indices)
    call test_report%assert(all(a_sorted == [-2._dp, -1._dp, 0._dp, 0._dp, 1._dp, 2._dp]), &
        message="Expect array ordered lowest to highest")

  end subroutine

  subroutine test_sort_index_2d(test_report)
    type(unit_test_type), intent(inout) :: test_report

    integer :: a(3,4), a_sorted(3, 4), expected_a_sorted(3, 4)
    integer :: i, j
    integer :: index_map(4)

    ! Each column corresponds to a multi-digital number i.e. [104, 113, 132, 121]
    a = transpose(reshape([1, 1, 1, 1, &
                           0, 1, 3, 2, &
                           4, 3, 2, 1 ], [4, 3]))


    index_map = sort_index_2d(size(a, 1), size(a, 2), a, size(a, 1))

    call test_report%assert(all(index_map ==[1, 2, 4, 3]))

    ! [104, 113, 132, 121] --> [104, 113, 121, 132]
    expected_a_sorted = transpose(reshape([1, 1, 1, 1, &
                                           0, 1, 2, 3, &
                                           4, 3, 1, 2 ], [4, 3]))

    do i = 1, 4
        j = index_map(i)
        a_sorted(:, j) = a(:, i)
    end do

    call test_report%assert(all(a_sorted == expected_a_sorted), &
      message="Expect `a` to be sorted according to the numbers defined by &
              concatenating each row, per column")

  end subroutine

end module
