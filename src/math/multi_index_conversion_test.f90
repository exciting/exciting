module multi_index_conversion_test
  use modmpi, only: mpiinfo
  use unit_test_framework, only : unit_test_type

  use math_utils, only: mod1
  use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices

  implicit none

  private
  public :: multi_index_conversion_test_driver
    
  
contains

  !> Run tests for multi index conversion
  subroutine multi_index_conversion_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure  
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 9

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    
    call test_multi_index_to_index(test_report)

    call test_index_to_multi_index(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('multi_index_conversion', kill_on_failure)
    else
      call test_report%report('multi_index_conversion')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine multi_index_conversion_test_driver


  !> Test index
  subroutine test_multi_index_to_index(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    integer :: i, j, k, l, m, composite_index, N(5), cnt
    integer, allocatable :: counters(:), composite_indices(:)

    N = [7, 2, 34, 23, 2]

    allocate(composite_indices(product(N)))
    allocate(counters(product(N)))

    ! loop over a multi index, going from the outer to inner corresponds
    ! to going from right to left.
    cnt = 0
    do m=1, N(5)
      do l=1, N(4)
        do k=1, N(3)
          do j=1, N(2)
            do i=1, N(1)
              cnt = cnt + 1
              composite_indices(cnt) = indices_to_composite_index([i, j, k, l, m], N)
              counters(cnt) = cnt
            end do 
          end do 
        end do 
      end  do 
    end do

    call test_report%assert(all(counters == composite_indices), &
                           ' Composite index returned disagrees with the loop counter cnt running from [1: product(N)]')
    end subroutine test_multi_index_to_index


  !> Test index
  subroutine test_index_to_multi_index(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    integer :: i, j, k, l, m, o, i_1, i_2, i_3, i_4, i_5, composite_index, cnt
    integer, allocatable :: multi_index(:), multi_index_test(:), N(:)
    logical, allocatable :: is_multi_index(:)

    ! Test index_to_double_index
    N = [7, 2]
    allocate(is_multi_index(product(N)))

    ! loop over a multi index, going from the outer to inner cooresponds
    ! to going from right to left.
    composite_index = 1
    do j=1, N(2)
      do i=1, N(1)
        call composite_index_to_indices(composite_index, N(1), i_1, i_2)
        is_multi_index(composite_index) = all([i, j] == [i_1, i_2])
        composite_index = composite_index + 1
      end do
    end do
    call test_report%assert(all(is_multi_index), &
                           'Multi indices returned by composite_index_to_indices does disagree with the nested indices i and j.')


    ! Test index_to_multi_index_shape
    N = [7, 2, 34, 23, 2, 8]
    allocate(multi_index_test(6))
    deallocate(is_multi_index)
    allocate(is_multi_index(product(N)))

    ! loop over a multi index, going from the outer to inner cooresponds
    ! to going from right to left.
    composite_index = 1
    do o=1, N(6)
      do m=1, N(5)
        do l=1, N(4)
          do k=1, N(3)
            do j=1, N(2)
              do i=1, N(1)
                multi_index = [i, j, k, l, m, o]
                call composite_index_to_indices(composite_index, N, multi_index_test)
                is_multi_index(composite_index) = all(multi_index_test == multi_index)
                composite_index = composite_index + 1
              end do 
            end do 
          end do 
        end  do 
      end do
    end do
    call test_report%assert(all(is_multi_index), &
                           'Multi indices returned by composite_index_to_indices does disagree with the nested indices i, j, k, l, m and o.')

    ! Test index_to_double_index
    N = [7, 2]
    deallocate(is_multi_index)
    allocate(is_multi_index(product(N)))

    ! loop over a multi index, going from the outer to inner cooresponds
    ! to going from right to left.
    composite_index = 1
    do j=1, N(2)
      do i=1, N(1)
        call composite_index_to_indices(composite_index, N, i_1, i_2)
        is_multi_index(composite_index) = all([i, j] == [i_1, i_2])
        composite_index = composite_index + 1
      end do
    end do
    call test_report%assert(all(is_multi_index), &
                            'Multi indices returned by composite_index_to_indices does disagree with the nested indices i and j.')

    ! Test index_to_triple_index
    N = [2, 3, 4]
    deallocate(is_multi_index)
    allocate(is_multi_index(product(N)))

    ! loop over a multi index, going from the outer to inner cooresponds
    ! to going from right to left.
    composite_index = 1
    do k=1, N(3)
      do j=1, N(2)
        do i=1, N(1)
          call composite_index_to_indices(composite_index, N, i_1, i_2, i_3)
          is_multi_index(composite_index) = all([i, j, k] == [i_1, i_2, i_3])
          composite_index = composite_index + 1
        end do
      end do
    end do
    call test_report%assert(all(is_multi_index), &
                            'Multi indices returned by composite_index_to_indices does disagree with the nested indices i, j and k.')


    ! Test index_to_quartuple_index
    N = [13, 21, 2, 5]
    deallocate(is_multi_index)
    allocate(is_multi_index(product(N)))

    ! loop over a multi index, going from the outer to inner cooresponds
    ! to going from right to left.
    composite_index = 1
    do l=1, N(4)
      do k=1, N(3)
        do j=1, N(2)
          do i=1, N(1)
            call composite_index_to_indices(composite_index, N, i_1, i_2, i_3, i_4)
            is_multi_index(composite_index) = all([i, j, k, l] == [i_1, i_2, i_3, i_4])
            composite_index = composite_index + 1
          end do
        end do
      end do
    end do
    call test_report%assert(all(is_multi_index), &
                            'Multi indices returned by composite_index_to_indices does disagree with the nested indices i, j, k and l.')


    ! Test index_to_quintuple_index
    N = [12, 5, 12, 3, 32]
    deallocate(is_multi_index)
    allocate(is_multi_index(product(N)))

    ! loop over a multi index, going from the outer to inner cooresponds
    ! to going from right to left.
    composite_index = 1
    do m=1, N(5)
      do l=1, N(4)
        do k=1, N(3)
          do j=1, N(2)
            do i=1, N(1)
              call composite_index_to_indices(composite_index, N, i_1, i_2, i_3, i_4, i_5)
              is_multi_index(composite_index) = all([i, j, k, l, m] == [i_1, i_2, i_3, i_4, i_5])
              composite_index = composite_index + 1
            end do
          end do
        end do
      end do
    end do
    call test_report%assert(all(is_multi_index), &
                            'Multi indices returned by composite_index_to_indices does disagree with the nested indices i, j, k, l and m.')
  end subroutine test_index_to_multi_index

end module multi_index_conversion_test



