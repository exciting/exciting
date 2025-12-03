! Created by  on 18.05.22.

module regular_grid_test
  use modmpi, only: mpiinfo
  use precision, only: dp
  use unit_test_framework, only : unit_test_type
  use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices
  use math_utils, only: all_close, fractional_part, identity_real_dp
  use grid_utils, only: grid_3d, mesh_1d, fft_frequencies
  use regular_grid, only: regular_grid_type, regular_grid_type
  use bravais_lattice, only: body_centered_orthorhombic, simple_monoclinic, simple_cubic

  implicit none

  private
  public :: regular_grid_test_driver


  contains

  !> Run tests for regular grid
  subroutine regular_grid_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Run and assert tests
    call test_regular_grid(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('regular_grid_test', kill_on_failure)
    else
      call test_report%report('regular_grid_test')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine regular_grid_test_driver


  !> Static tests for regular grid
  subroutine test_regular_grid(test_report)
    !> Unit test object
    type(unit_test_type), intent(inout) :: test_report

    ! 3 x 2 x 4 mesh with [1, 1, 1] as first point
    integer :: ref_mesh_0(3, 24) =  reshape(&
      [ 1, 1, 1,  2, 1, 1,  3, 1, 1, &
        1, 2, 1,  2, 2, 1,  3, 2, 1, &

        1, 1, 2,  2, 1, 2,  3, 1, 2, &
        1, 2, 2,  2, 2, 2,  3, 2, 2, &

        1, 1, 3,  2, 1, 3,  3, 1, 3, &
        1, 2, 3,  2, 2, 3,  3, 2, 3, &

        1, 1, 4,  2, 1, 4,  3, 1, 4, &
        1, 2, 4,  2, 2, 4,  3, 2, 4], [3, 24])

    ! 3 x 2 x 4 mesh with [3, 1, 2] as first point
    integer, parameter :: ref_mesh_1(3, 24) = reshape( &
      [ 3, 1, 2,  1, 1, 2,  2, 1, 2, &
        3, 2, 2,  1, 2, 2,  2, 2, 2, &

        3, 1, 3,  1, 1, 3,  2, 1, 3, &
        3, 2, 3,  1, 2, 3,  2, 2, 3, &

        3, 1, 4,  1, 1, 4,  2, 1, 4, &
        3, 2, 4,  1, 2, 4,  2, 2, 4, &

        3, 1, 1,  1, 1, 1,  2, 1, 1, &
        3, 2, 1,  1, 2, 1,  2, 2, 1], [3, 24])

    ! 3 x 2 x 4 grid, starting at integer coordinate [1, 1, 1], with an offset of [0.24, 0.03, 0.38]
    real(dp), parameter :: ref_grid_0(3, 24) = reshape( &
      [ 0.24_dp, 0.03_dp, 0.38_dp,   1.24_dp, 0.03_dp, 0.38_dp,   2.24_dp, 0.03_dp, 0.38_dp, &
        0.24_dp, 1.03_dp, 0.38_dp,   1.24_dp, 1.03_dp, 0.38_dp,   2.24_dp, 1.03_dp, 0.38_dp, &

        0.24_dp, 0.03_dp, 1.38_dp,   1.24_dp, 0.03_dp, 1.38_dp,   2.24_dp, 0.03_dp, 1.38_dp, &
        0.24_dp, 1.03_dp, 1.38_dp,   1.24_dp, 1.03_dp, 1.38_dp,   2.24_dp, 1.03_dp, 1.38_dp, &

        0.24_dp, 0.03_dp, 2.38_dp,   1.24_dp, 0.03_dp, 2.38_dp,   2.24_dp, 0.03_dp, 2.38_dp, &
        0.24_dp, 1.03_dp, 2.38_dp,   1.24_dp, 1.03_dp, 2.38_dp,   2.24_dp, 1.03_dp, 2.38_dp, &

        0.24_dp, 0.03_dp, 3.38_dp,   1.24_dp, 0.03_dp, 3.38_dp,   2.24_dp, 0.03_dp, 3.38_dp, &
        0.24_dp, 1.03_dp, 3.38_dp,   1.24_dp, 1.03_dp, 3.38_dp,   2.24_dp, 1.03_dp, 3.38_dp ], [3, 24])


    ! 3 x 2 x 4 grid, starting at integer coordinate [3, 1, 2], with an offset of [0.24, 0.03, 0.38]
    real(dp), parameter :: ref_grid_1(3, 24) = reshape( &
      [ 2.24_dp, 0.03_dp, 1.38_dp,   0.24_dp, 0.03_dp, 1.38_dp,   1.24_dp, 0.03_dp, 1.38_dp, &
        2.24_dp, 1.03_dp, 1.38_dp,   0.24_dp, 1.03_dp, 1.38_dp,   1.24_dp, 1.03_dp, 1.38_dp, &

        2.24_dp, 0.03_dp, 2.38_dp,   0.24_dp, 0.03_dp, 2.38_dp,   1.24_dp, 0.03_dp, 2.38_dp, &
        2.24_dp, 1.03_dp, 2.38_dp,   0.24_dp, 1.03_dp, 2.38_dp,   1.24_dp, 1.03_dp, 2.38_dp, &

        2.24_dp, 0.03_dp, 3.38_dp,   0.24_dp, 0.03_dp, 3.38_dp,   1.24_dp, 0.03_dp, 3.38_dp, &
        2.24_dp, 1.03_dp, 3.38_dp,   0.24_dp, 1.03_dp, 3.38_dp,   1.24_dp, 1.03_dp, 3.38_dp, &

        2.24_dp, 0.03_dp, 0.38_dp,   0.24_dp, 0.03_dp, 0.38_dp,   1.24_dp, 0.03_dp, 0.38_dp, &
        2.24_dp, 1.03_dp, 0.38_dp,   0.24_dp, 1.03_dp, 0.38_dp,   1.24_dp, 1.03_dp, 0.38_dp ], [3, 24])

    ! Order of the grid points (ascending from 1 to 24 with a spacing of 1)
    integer, parameter :: ref_order(24) = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, &
                                           13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

    type(regular_grid_type) :: grid

    integer :: sampling(3), multi_index_first(3)
    real(dp) :: offset(3), volume_element(3, 3)
    real(dp), allocatable :: sub_grid(:, :)

    sampling = [3, 2, 4]
    multi_index_first = [1, 1, 1]
    offset = [.24_dp, .03_dp, .38_dp]
    volume_element = reshape( [1.0_dp, 0.0_dp, 0.0_dp, &
                               0.0_dp, 1.0_dp, 0.0_dp, &
                               0.0_dp, 0.0_dp, 1.0_dp], [3, 3])

    grid = regular_grid_type(sampling, multi_index_first, offset, volume_element)

    ! composite_index -> coord
    call test_report%assert( &
      all_close(grid%coordinate(9), [2.24_dp, 0.03_dp, 1.38_dp]), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%coordinate_array(9) /= [2.24_dp, 0.03_dp, 1.38_dp].' &
      )
    call test_report%assert( &
      all_close(grid%coordinate_array(ref_order), ref_grid_0), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%coordinate_array(ref_order)/= ref_grid.' &
      )

    ! multi_index -> coord
    call test_report%assert( &
      all_close(grid%coordinate([3, 2, 1]), [2.24_dp, 1.03_dp, 0.38_dp]), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%coordinate([3, 2, 1])/= [2.24_dp, 1.03_dp, 0.38_dp].' &
      )
    call test_report%assert( &
      all_close(grid%coordinate_array(ref_mesh_0), ref_grid_0), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%coordinate_array(ref_mesh_0) /= ref_grid_0.' &
      )

    call test_report%assert( &
        all_close(grid%coordinate_array(), ref_grid_0), &
        'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
            grid%coordinate_array() is not the same as ref_grid_0' &
        )

    ! coord -> multi_index
    call test_report%assert( &
      all(grid%multi_index([1.24_dp, 1.03_dp, 1.38_dp]) == [2, 2, 2]), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%multi_index([1.24_dp, 1.03_dp, 1.38_dp]) /= [2, 2, 2].' &
      )
    call test_report%assert( &
      all(grid%multi_index(ref_grid_0) == ref_mesh_0), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%multi_index(ref_grid_0) /= ref_mesh_0.' &
      )

    ! coord -> composite_index
    call test_report%assert( &
      grid%composite_index([2.24_dp, 0.03_dp, 1.38_dp]) == 9, &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%composite_index(2.24_dp, 0.03_dp, 1.38_dp) /= 9.' &
      )
    call test_report%assert( &
      all(grid%composite_index(ref_grid_0) == ref_order), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        grid%composite_index(ref_grid_0) /= ref_order.' &
      )

    ! test is_on_grid
    call test_report%assert(grid%is_on_grid([1.24_dp, 1.03_dp, 1.38_dp]), &
        'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
            this%is_on_grid([1.24_dp, 1.03_dp, 1.38_dp]) is false.')

    call test_report%assert(.not. grid%is_on_grid([1.56_dp, 1.03_dp, 1.38_dp]), &
        'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
            this%is_on_grid([1.56_dp, 1.03_dp, 1.38_dp]) is true.')

    call test_report%assert(all(grid%is_on_grid(reshape([1.24_dp, 1.03_dp, 1.38_dp, &
                                                         1.56_dp, 1.03_dp, 1.38_dp], [3, 2])) .eqv. [.true., .false.]), &
        'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
            this%is_on_grid(reshape([1.24_dp, 1.03_dp, 1.38_dp, &
                                     1.56_dp, 1.03_dp, 1.38_dp], [3, 2])) is [.false., .true.]')


    sampling = [3, 2, 4]
    multi_index_first = [3, 1, 2]
    offset = [.24_dp, .03_dp, .38_dp]
    volume_element = reshape( [1.0_dp, 0.0_dp, 0.0_dp, &
                               0.0_dp, 1.0_dp, 0.0_dp, &
                               0.0_dp, 0.0_dp, 1.0_dp ], [3, 3])

    grid = regular_grid_type(sampling, multi_index_first, offset, volume_element)

    ! composite_index -> coord
    call test_report%assert( &
      all_close(grid%coordinate(9), [1.24_dp, 0.03_dp, 2.38_dp]), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%coordinate_array(9) /= [1.24_dp, 0.03_dp, 2.38_dp].' &
      )
    call test_report%assert( &
      all_close(grid%coordinate_array(ref_order), ref_grid_1), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%coordinate_array(ref_order)/= ref_grid_1.' &
      )

    ! multi_index -> coord
    call test_report%assert( &
      all_close(grid%coordinate([3, 2, 1]), [2.24_dp, 1.03_dp, 0.38_dp]), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%coordinate([3, 2, 1])/= [2.24_dp, 1.03_dp, 0.38_dp].' &
      )
    call test_report%assert( &
      all_close(grid%coordinate_array(ref_mesh_1), ref_grid_1), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%coordinate_array(ref_mesh_1) /= ref_grid_1.' &
      )

    call test_report%assert( &
        all_close(grid%coordinate_array(), ref_grid_1), &
        'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
            grid%coordinate_array() is not the same as ref_grid_1.' &
        )

    ! coord -> multi_index
    call test_report%assert( &
      all(grid%multi_index([1.24_dp, 1.03_dp, 1.38_dp]) == [2, 2, 2]), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%multi_index([1.24_dp, 1.03_dp, 1.38_dp]) /= [2, 2, 2].' &
      )
    call test_report%assert( &
      all(grid%multi_index(ref_grid_1) == ref_mesh_1), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%multi_index(ref_grid_1) /= ref_mesh_1.' &
      )

    ! coord -> composite_index
    call test_report%assert( &
      grid%composite_index([2.24_dp, 0.03_dp, 1.38_dp]) == 1, &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%composite_index([2.24_dp, 0.03_dp, 1.38_dp]) /= 1.' &
      )
    call test_report%assert( &
      all(grid%composite_index(ref_grid_1) == ref_order), &
      'regular_grid static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        grid%composite_index(ref_grid_1) /= ref_order.' &
      )
  end subroutine test_regular_grid

end module regular_grid_test