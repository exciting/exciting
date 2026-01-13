! Created by Benedikt Maurer on 31.07.22.

module regular_mesh_test
  use modmpi, only: mpiinfo
  use precision, only: str_512
  use unit_test_framework, only : unit_test_type
  use math_utils, only: mod1
  use grid_utils, only: mesh_1d
  use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices
  use regular_mesh, only: regular_mesh_type

  implicit none

  private
  public :: regular_mesh_test_driver

contains

  !> Run tests for regular_mesh
  subroutine regular_mesh_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report

    ! Initialize test object
    call test_report%init( mpiglobal)

    ! Static test for 3 x 2 x 4 mesh
    call test_regular_grid_static(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('regular_mesh_test', kill_on_failure)
    else
      call test_report%report('regular_mesh_test')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine regular_mesh_test_driver

  subroutine test_regular_grid_static(test_report)
    !> Unit test object
    type(unit_test_type), intent(inout) :: test_report

    ! 3 x 2 x 4 mesh with [1, 1, 1] as first point
    integer :: ref_mesh_0(3, 24) =  reshape( &
      [ 1, 1, 1,  2, 1, 1,  3, 1, 1, &
        1, 2, 1,  2, 2, 1,  3, 2, 1, &

        1, 1, 2,  2, 1, 2,  3, 1, 2, &
        1, 2, 2,  2, 2, 2,  3, 2, 2, &

        1, 1, 3,  2, 1, 3,  3, 1, 3, &
        1, 2, 3,  2, 2, 3,  3, 2, 3, &

        1, 1, 4,  2, 1, 4,  3, 1, 4, &
        1, 2, 4,  2, 2, 4,  3, 2, 4 ], [3, 24])

    ! 3 x 2 x 4 mesh with [3, 1, 2] as first point
    integer, parameter :: ref_mesh_1(3, 24) = reshape( &
      [ 3, 1, 2,  1, 1, 2,  2, 1, 2, &
        3, 2, 2,  1, 2, 2,  2, 2, 2, &

        3, 1, 3,  1, 1, 3,  2, 1, 3, &
        3, 2, 3,  1, 2, 3,  2, 2, 3, &

        3, 1, 4,  1, 1, 4,  2, 1, 4, &
        3, 2, 4,  1, 2, 4,  2, 2, 4, &

        3, 1, 1,  1, 1, 1,  2, 1, 1, &
        3, 2, 1,  1, 2, 1,  2, 2, 1 ], [3, 24])

    ! Order of the mesh points (ascending from 1 to 24 with a spacing of 1)
    integer, parameter :: ref_order(24) = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, &
                                           13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

    ! Alternative mapping
    integer, parameter :: index_map(8) = [ 1, 4, 6, 2, 20, 8, 14, 13]

    type(regular_mesh_type) :: mesh

    ! Test 3 x 2 x 4 mesh with [1, 1, 1] as first integer position
    mesh = regular_mesh_type([3, 2, 4], [1, 1, 1])

    ! Test mesh%number_of_points
    call test_report%assert( &
      mesh%number_of_points() == 24, &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%number_of_points() /= 1." &
      )

    ! Test mesh%multi_index
    call test_report%assert( &
      all(mesh%multi_index(5) == [2, 2, 1]), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%multi_index(5) /= [2, 2, 1])" &
      )
    call test_report%assert( &
      all(mesh%multi_index(ref_order) == ref_mesh_0), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%multi_index(ref_order) /= ref_mesh_0." &
      )
    call test_report%assert( &
      all(mesh%multi_index() == ref_mesh_0), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%multi_index() /= ref_mesh_0." &
      )

    ! Test mesh%composite_index
    call test_report%assert( &
       mesh%composite_index([2, 2, 3]) == 17, &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%composite_index([2, 2, 3]) /= 17." &
      )
    call test_report%assert( &
      all(mesh%composite_index(ref_mesh_0) == ref_order), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%composite_index(ref_mesh_0) /= ref_order." &
      )
    call test_report%assert( &
      all(mesh%composite_index() == ref_order), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [1, 1, 1]: &
        mesh%composite_index() /= ref_order." &
      )

    ! Test 3 x 2 x 4 mesh with [1, 1, 1] as first integer position
    mesh = regular_mesh_type([3, 2, 4], [3, 1, 2])

    ! Test mesh%number_of_points
    call test_report%assert( &
      mesh%number_of_points() == 24, &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%number_of_points() /= 1." &
      )

    ! Test mesh%multi_index
    call test_report%assert( &
      all(mesh%multi_index(5) == [1, 2, 2]), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%multi_index(5) /= [1, 2, 2])" &
      )
    call test_report%assert( &
      all(mesh%multi_index(ref_order) == ref_mesh_1), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%multi_index(ref_order) /= ref_mesh." &
      )
    call test_report%assert( &
      all(mesh%multi_index() == ref_mesh_1), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%multi_index() /= ref_mesh_1." &
      )

    ! Test mesh%composite_index
    call test_report%assert( &
      mesh%composite_index([2, 2, 3]) == 12, &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%composite_index([2, 2, 3]) /= 17." &
      )
    call test_report%assert( &
      all(mesh%composite_index(ref_mesh_1) == ref_order), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%composite_index(ref_mesh) /= ref_order." &
      )
    call test_report%assert( &
      all(mesh%composite_index() == ref_order), &
      "regular_mesh static test: sampling = [3, 2, 4], multi_index_first = [3, 1, 2]: &
        mesh%composite_index() /= ref_order." &
      )


  end subroutine

end module regular_mesh_test