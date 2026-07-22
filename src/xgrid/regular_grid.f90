! Created by  on 09/02/2022.
!> Module for the base class of a regular grid. Any implementation of a new regular grid type can extend this
!> class.
module regular_grid
#include "asserts.fpp"
  use precision, only: dp
  use math_utils, only: mod1, all_close, fractional_part, integer_part
  use grid_utils, only: mesh_1d
  use linear_algebra_3d, only: determinant_3d, inverse_3d
  use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices
  use regular_mesh, only: regular_mesh_type

  implicit none

  private
  public :: regular_grid_type

  !> Default tolerance to decide if a double is zero
  real, parameter :: default_tolerance = 1e-10_dp

  !> Class for managing regular grids. Allows for switching between the index `composite_index`, the integer coordinate `multi_index` and
  !> the point coordinate `coordinate` easily. The class is extends [[regular_mesh/regular_mesh_type]]. For details of the
  !> `composite_index <-> multi_index` conversion refer to the [[regular_mesh]] documentation.
  !>
  !> The transformation from integer coordinates to point coordinates is defined as
  !> \[
  !>   \mathbf{r} = \mathbf{V} \cdot \left( \mathbf{i_{xyz} - 1} + \mathbf{o} \right),
  !> \]
  !> where \( \mathbf{V} \) is a paralleliped volume element, defined column-wise by it's span vectors,
  !> \(\mathbf{i_{xyz}\) is the integer position of the point and \( \mathbf{o} \) the offset of the grid.
  !> \( \mathbf{V} \) and \( \mathbf{o} \) are global for the grid.
  type, extends(regular_mesh_type) :: regular_grid_type
    !> Offset of the real coordinates
    real(dp) :: offset(3)
    !> Volume element of the grid
    real(dp) :: volume_element(3, 3)
    !> Tolerance to decide if a double is zero
    real(dp) :: tol = default_tolerance

    contains

    procedure :: is_on_grid_single, is_on_grid_array
    generic :: is_on_grid => is_on_grid_single, is_on_grid_array

    procedure :: coordinate_to_composite_index, coordinate_to_composite_index_array
    generic :: composite_index => coordinate_to_composite_index, coordinate_to_composite_index_array

    procedure :: coordinate_to_multi_index, coordinate_to_multi_indices
    generic :: multi_index => coordinate_to_multi_index, coordinate_to_multi_indices

    procedure :: composite_index_to_coordinate, multi_index_to_coordinate
    generic :: coordinate => composite_index_to_coordinate, multi_index_to_coordinate

    procedure :: composite_index_to_coordinate_array, multi_index_to_coordinate_array, full_coordinate_array
    generic :: coordinate_array => multi_index_to_coordinate_array, composite_index_to_coordinate_array, full_coordinate_array
  end type regular_grid_type

  interface regular_grid_type
    procedure :: setup_regular_grid_type
  end interface regular_grid_type

  contains

  !> Default constructor
  function setup_regular_grid_type(sampling, multi_index_first, offset, volume_element, tol) result(this)
    !> Number of grid points per dimension.
    integer, intent(in) :: sampling(3)
    !> Integer coordinate of the first point.
    integer, intent(in) :: multi_index_first(3)
    !> Offset of the coordinates, in lattice coordinates
    real(dp), intent(in) :: offset(3)
    !> Volume element of the grid
    real(dp), intent(in) :: volume_element(3, 3)
    !> Tolerance to decide if a double is zero
    real(dp), intent(in), optional :: tol

    type(regular_grid_type) :: this

    real(dp) :: tol_local

    tol_local = default_tolerance
    if (present(tol)) tol_local = tol

    CALL_ASSERT(all(sampling > 0), 'sampling <= 0.')
    CALL_ASSERT(all(multi_index_first > 0), 'multi_index_first <= 0.')
    CALL_ASSERT(all(multi_index_first <= sampling), 'multi_index_first > sampling.')
    CALL_ASSERT(tol_local >= 0._dp, 'tol < 0._dp.')

    this%sampling =  sampling
    this%multi_index_first = multi_index_first
    this%offset = offset
    this%volume_element = volume_element
    this%tol = tol_local
  end function setup_regular_grid_type

  !> Verify that a coordinate is on the grid.
  logical function is_on_grid_single(this, coordinate)
    !> Regular grid instance.
    class(regular_grid_type), intent(in) :: this
    !> Coordinate to check on.
    real(dp), intent(in) :: coordinate(3)

    integer :: multi_index(3)
    real(dp) :: coordinate_local(3)

    coordinate_local = matmul(inverse_3d(this%volume_element), coordinate) - this%offset
    multi_index = idnint(coordinate_local) + 1
    is_on_grid_single = all_close(this%coordinate(multi_index), coordinate, this%tol)
  end function is_on_grid_single

  !> Verify that a set of coordinates is on the grid. The output is a logical array with the corresponding result.
  function is_on_grid_array(this, coordinate_array) result(is_on_grid_results)
    !> Regular grid instance.
    class(regular_grid_type), intent(in) :: this
    !> Coordinates to check on.
    real(dp), intent(in) :: coordinate_array(:, :)

    logical, allocatable :: is_on_grid_results(:)

    integer :: n_coordinates, i

    CALL_ASSERT(size(coordinate_array, 1) == 3, 'size(coordinate_array, 1) /= 3.')

    n_coordinates = size(coordinate_array, 2)

    allocate(is_on_grid_results(n_coordinates))
    do i=1, n_coordinates
      is_on_grid_results(i) = is_on_grid_single(this, coordinate_array(:, i))
    end do
  end function is_on_grid_array

  !---------------------------------------------------------------------------------------------
  ! Coordinate access
  !---------------------------------------------------------------------------------------------

  ! Multi index to coordinate

  !> Transform integer position to point coordinate.
  function multi_index_to_coordinate(this, multi_index) result(coordinate)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Inter position
    integer, intent(in) :: multi_index(3)

    real(dp) :: coordinate(3)

    coordinate = matmul(this%volume_element, real(multi_index - 1, dp) + this%offset)
  end function multi_index_to_coordinate

  !> Transform array of integer positions to array point coordinates.
  function multi_index_to_coordinate_array(this, multi_indices) result(coordinate_array)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Array of integer position (column-wise)
    integer, intent(in) :: multi_indices(:, :)

    real(dp), allocatable :: coordinate_array(:, :)

    integer :: number_of_points, i

    CALL_ASSERT(size(multi_indices, 1) == 3, 'size(multi_indices, 1) /= 3')

    number_of_points = size(multi_indices, 2)
    allocate(coordinate_array(3, number_of_points))

    do i=1, number_of_points
      coordinate_array(:, i) = multi_index_to_coordinate(this, multi_indices(:, i))
    end do
  end function multi_index_to_coordinate_array

  ! Composite index to coordinate

  !> Transform point index to point coordinate
  function composite_index_to_coordinate(this, composite_index) result(coordinate)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Point index
    integer, intent(in) :: composite_index

    real(dp) :: coordinate(3)

    integer :: composite_index_local, multi_index(3)

    ! Save some routine calls:
    ! Same implementation as regular_mesh/composite_index_to_multi_index to obtain multi_index
    composite_index_local = mod1(composite_index, this%number_of_points())
    call composite_index_to_indices(composite_index_local, this%sampling, multi_index)
    multi_index = mod1(multi_index + this%multi_index_first - 1, this%sampling)

    ! Same implementation as multi_index_to_coordinate to obatain coordinate
    coordinate = matmul(this%volume_element, real(multi_index - 1, dp) + this%offset)
  end function composite_index_to_coordinate

  !> Transform array of point indices to array of coordinates
  function composite_index_to_coordinate_array(this, composite_indices) result(coordinate_array)
    !> Regular grid instance.
    class(regular_grid_type), intent(in) :: this
    !> Array of point indices.
    integer, intent(in) :: composite_indices(:)

    real(dp), allocatable :: coordinate_array(:, :)

    integer :: number_of_points, i

    number_of_points = size(composite_indices)
    allocate(coordinate_array(3, number_of_points))

    do i=1, number_of_points
      coordinate_array(:, i) = composite_index_to_coordinate(this, composite_indices(i))
    end do
  end function composite_index_to_coordinate_array

  ! Access all coordinates

  !> Return all coordinates of the grid.
  function full_coordinate_array(this) result(coordinate_array)
    !> Regular grid instance.
    class(regular_grid_type), intent(in) :: this

    real(dp), allocatable :: coordinate_array(:, :)

    coordinate_array = this%coordinate_array(mesh_1d(1, this%number_of_points()))
  end function full_coordinate_array

  !---------------------------------------------------------------------------------------------
  ! Multi index access
  !---------------------------------------------------------------------------------------------

  !> Transform point coordinate to integer position
  function coordinate_to_multi_index(this, coordinate) result(multi_index)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Point coordinate
    real(dp), intent(in) :: coordinate(3)

    integer :: multi_index(3)

    real(dp) :: coordinate_local(3)

    CALL_ASSERT(this%is_on_grid(coordinate), 'coordinate is not on grid.')

    coordinate_local = matmul(inverse_3d(this%volume_element), coordinate) - this%offset
    multi_index = idnint(coordinate_local) + 1
  end function coordinate_to_multi_index

  !> Transform array of point coordinates to array of integer positions.
  function coordinate_to_multi_indices(this, coordinate_array) result(multi_indices)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Array of point coordinate (column-wise)
    real(dp), intent(in) :: coordinate_array(:, :)

    integer, allocatable :: multi_indices(:, :)

    integer :: number_of_points, i

    CALL_ASSERT(size(coordinate_array, 1) == 3, 'size(coordinate_array, 1) /= 3.')

    number_of_points = size(coordinate_array, 2)
    allocate(multi_indices(3, number_of_points))

    do i=1, number_of_points
      multi_indices(:, i) = coordinate_to_multi_index(this, coordinate_array(:, i))
    end do
  end function coordinate_to_multi_indices

  !---------------------------------------------------------------------------------------------
  ! Composite index access
  !---------------------------------------------------------------------------------------------

  !> Transform point coordinate to point index
  function coordinate_to_composite_index(this, coordinate) result(composite_index)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Point coordinate
    real(dp), intent(in) :: coordinate(3)

    integer :: composite_index

    integer :: multi_index(3)
    real(dp) :: coordinate_local(3)

    !CALL_ASSERT(this%is_on_grid(coordinate), 'coordinate is not on grid.')
    if (.not. this%is_on_grid(coordinate)) then
      composite_index = -1
      return 
    end if  

    ! Same implementation as coordinate_to_multi_index to obtain multi_index
    coordinate_local = matmul(inverse_3d(this%volume_element), coordinate) - this%offset
    multi_index = idnint(coordinate_local) + 1

    ! Same implementation as multi_index_to_composite_index to obatain composite_index
    multi_index = mod1(multi_index - this%multi_index_first + 1, this%sampling)
    composite_index = indices_to_composite_index(multi_index, this%sampling)
  end function coordinate_to_composite_index

  !> Transform array of point coordinates to array of point indices
  function coordinate_to_composite_index_array(this, coordinate_array) result(composite_indices)
    !> Regular grid instance
    class(regular_grid_type), intent(in) :: this
    !> Array  of point coordinates
    real(dp), intent(in) :: coordinate_array(:, :)

    integer, allocatable :: composite_indices(:)

    integer :: number_of_points, i

    CALL_ASSERT(size(coordinate_array, 1) == 3, 'size(coordinate_array, 1) /= 3.')

    number_of_points = size(coordinate_array, 2)
    allocate(composite_indices(number_of_points))

    do i=1, number_of_points
      composite_indices(i) = coordinate_to_composite_index(this, coordinate_array(:, i))
    end do
  end function coordinate_to_composite_index_array

end module regular_grid