! Created by  on 17/02/2022.

module regular_grid_constructors_test
  use modmpi, only: mpiinfo
  use precision, only: dp
  use unit_test_framework, only : unit_test_type
  use regular_grid, only: regular_grid_type
  use  xgrid, only: setup_unitcell_grid, setup_fft_grid
  use bravais_lattice, only: rhombohedral_hex_setting
  use grid_utils, only: mesh_1d, grid_3d, fft_frequencies
  use math_utils, only: all_close
  use unit_cell_utils, only: reciprocal_lattice

  implicit none

  private
  public :: regular_grid_constructors_test_driver

  contains

  !> Run tests for regular grid
  subroutine regular_grid_constructors_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails
    logical, optional :: kill_on_failure
    !> test object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 3

    ! Initialize test object
    call test_report%init(n_assertions, mpiglobal)

    ! Run and assert tests
    call test_setup_unitcell_grid(test_report)

    call test_setup_fft_grid(test_report)

    ! report results
    if (present(kill_on_failure)) then
      call test_report%report('regular_G_grid', kill_on_failure)
    else
      call test_report%report('regular_G_grid')
    end if

    ! Finalise test object
    call test_report%finalise()
  end subroutine regular_grid_constructors_test_driver

  subroutine test_setup_unitcell_grid(test_report)
    type(unit_test_type) :: test_report

    type(regular_grid_type) :: unitcell_grid, test_fft_grid_1, test_fft_grid_2

    logical :: is_ascending
    integer :: sampling(3)
    integer, allocatable :: ref_order(:)
    real(dp) :: offset(3), lattice(3, 3)
    real(dp), allocatable :: ref_coords(:, :)

    sampling = [3, 4, 2]
    offset = [0._dp, 0._dp, 0._dp]
    lattice = rhombohedral_hex_setting(3._dp, 4._dp)

    ref_coords = matmul(lattice, grid_3d(sampling, include_last = .false.))
    ref_order = mesh_1d(1, product(sampling))

    unitcell_grid = setup_unitcell_grid(sampling, offset, lattice)

    call test_report%assert(all_close(unitcell_grid%coordinate_array(ref_order), ref_coords), &
        'unitcell_grid%coordinate_array(ref_order) is not ref_coords.')
  end subroutine test_setup_unitcell_grid

  subroutine test_setup_fft_grid(test_report)
    type(unit_test_type) :: test_report

    type(regular_grid_type) :: unitcell_grid, test_fft_grid_1, test_fft_grid_2

    logical :: is_ascending
    integer :: sampling(3)
    integer, allocatable :: ref_order(:)
    real(dp) :: radius, offset(3), lattice(3, 3)
    real(dp), allocatable :: ref_coords(:, :)

    sampling = [3, 4, 2]
    offset = [0._dp, 0._dp, 0._dp]
    lattice = rhombohedral_hex_setting(3._dp, 4._dp)

    ref_coords = matmul(reciprocal_lattice(lattice), fft_frequencies(sampling))
    ref_order = mesh_1d(1, product(sampling))

    unitcell_grid = setup_unitcell_grid(sampling, offset, lattice)
    test_fft_grid_1 = setup_fft_grid(sampling, offset, reciprocal_lattice(lattice))
    test_fft_grid_2 = setup_fft_grid(unitcell_grid)

    call test_report%assert(all_close(test_fft_grid_1%coordinate_array(ref_order), ref_coords), &
        'test_fft_grid_1%coordinate_array(ref_order) is not ref_coords.')

    call test_report%assert(all_close(test_fft_grid_2%coordinate_array(ref_order), ref_coords), &
        'test_fft_grid_2%coordinate_array(ref_order) is not ref_coords.')
  end subroutine test_setup_fft_grid

end module regular_grid_constructors_test