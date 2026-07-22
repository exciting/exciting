! Created by  on 21/09/2022.

module regular_grid_constructors
  use precision, only: dp
#include "asserts.fpp"
  use math_utils, only: all_close
  use grid_utils, only: mesh_1d
  use unit_cell_utils, only: reciprocal_lattice
  use regular_grid, only: regular_grid_type

  implicit none

  private
  public :: setup_unitcell_grid, setup_fft_grid

  !> Default tolerance to decide if a double is zero
  real, parameter :: default_tolerance = 1e-10_dp

  interface setup_fft_grid
    procedure :: setup_fft_frequency_grid_from_parameters, setup_fft_frequency_grid_from_unitcell_grid
  end interface setup_fft_grid

  contains

  !> Setup a regular grid that samples a unit cell.
  function setup_unitcell_grid(sampling, offset, lattice_vectors, tol) result(this)
    !> Number of grid points per dimension.
    integer, intent(in) :: sampling(3)
    !> Offset of the real coordinates
    real(dp), intent(in) :: offset(3)
    !> Column-wise lattice vectors of the unit cell.
    real(dp), intent(in) :: lattice_vectors(3, 3)
    !> Tolerance to decide if a double is zero
    real(dp), intent(in), optional :: tol

    type(regular_grid_type) :: this

    real(dp) :: tol_local

    tol_local = default_tolerance
    if (present(tol)) tol_local = tol

    CALL_ASSERT(all(offset >= 0._dp), 'offset < 0.0')
    CALL_ASSERT(all(offset < 1._dp / real(sampling, dp)), 'offset >= 1 / sampling.')

    this = regular_grid_type( &
              sampling = sampling, &
              multi_index_first = [1, 1, 1], &
              offset = offset, &
              volume_element = lattice_vectors / transpose(spread(real(sampling, dp), 2, 3)), &
              tol = tol_local &
        )
  end function setup_unitcell_grid

  !> Setup FFT frequency grid from a grid that samples a unit cell.
  function setup_fft_frequency_grid_from_parameters(sampling, offset, reciprocal_lattice_vectors, tol) result(this)
    !> Number of grid points per dimension.
    integer, intent(in) :: sampling(3)
    !> Offset of the real coordinates
    real(dp), intent(in) :: offset(3)
    !> Column-wise reciprocal lattice vectors of the unit cell.
    real(dp), intent(in) :: reciprocal_lattice_vectors(3, 3)
    !> Tolerance to decide if a double is zero
    real(dp), intent(in), optional :: tol

    type(regular_grid_type) :: this

    real(dp) :: tol_local

    tol_local = default_tolerance
    if (present(tol)) tol_local = tol

    CALL_ASSERT(all(offset >= 0._dp), 'offset < 0.0')
    CALL_ASSERT(all(offset < 1._dp ), 'offset >= 1.0')

    this = regular_grid_type( &
              sampling = sampling, &
              multi_index_first = floor(sampling / 2._dp) + 1, &
              offset = offset - real(floor(sampling / 2._dp), dp), &
              volume_element = reciprocal_lattice_vectors, &
              tol = tol_local &
        )
  end function setup_fft_frequency_grid_from_parameters

  function setup_fft_frequency_grid_from_unitcell_grid(unitcell_grid) result(this)
    !> Unit cell grid.
    type(regular_grid_type), intent(in) :: unitcell_grid

    type(regular_grid_type) :: this

    this = regular_grid_type( &
              sampling = unitcell_grid%sampling, &
              multi_index_first = floor(unitcell_grid%sampling / 2._dp) + 1, &
              offset = -real(floor(unitcell_grid%sampling / 2._dp), dp), &
              volume_element = reciprocal_lattice(unitcell_grid%volume_element &
                  * transpose(spread(unitcell_grid%sampling, 2, 3))), &
              tol = unitcell_grid%tol &
        )
  end function setup_fft_frequency_grid_from_unitcell_grid

end module regular_grid_constructors