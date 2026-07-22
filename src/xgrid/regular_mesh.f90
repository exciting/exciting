!> Module for generating regular meshes.
module regular_mesh
#include "asserts.fpp"
    use math_utils, only: mod1
    use grid_utils, only: mesh_1d
    use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices

    implicit none

    private
    public :: regular_mesh_type

    !> Regular mesh defined by its `sampling`, (number of grid points per dimension), and `multi_index_first`
    !> (multi index corresponding to `composite_index = 1`).
    !>
    !> Provides the mapping between `composite_index` and the `multi_index = [index_1, index_2, index_3]`. The ordering
    !> is column-major, thus the `index_1` changes the fastest with `index`. For details on the transformation, refer to
    !> [[multi_index_conversion/indices_to_composite_index]] and [[multi_index_conversion/composite_index_to_indices]].
    !>
    !> For example, let `sampling = [3, 2, 4]` and `multi_index_first = [1, 1, 1]`:
    !>
    !> - `1 --> [1, 1, 1]`
    !> - `2 --> [2, 1, 1]`
    !>
    !> `multi_index` is returned in the limits between `[1, 1, 1]' and `sampling`. The next point in a certain
    !> direction is obtained by adding `1` to the respective integer coordinate. If it is greater than the
    !> number of mesh points in that dimension, the counting restarts at 1 and the next integer coordinate is increased
    !> by 1. The relation 'next integer coordinate' is cyclic and follows
    !> `index_1 -> index_2 -> index_3 -> index_1 - index_2 -> ...`.
    !>
    !> For example, let `sampling = [3, 2, 4]` and `multi_index_first = [1, 1, 1]`:
    !>
    !> - Next neighbor in first direction: `[1, 2, 4] --> [2, 2, 4]`
    !> - Next neighbor in second direction: `[2, 2, 4] --> [2, 1, 1]`.
    !>
    !> `multi_index_first` allows for assigning any `multi_index' of the mesh to `index = 1`. The counting follows as
    !> described above.
    type regular_mesh_type
      !> Number of grid points per dimension.
      integer :: sampling(3)
      !> Multi index, assigned to `composite_index = 1` for indentity mapping.
      integer :: multi_index_first(3)

      contains

      procedure :: number_of_points

      procedure :: composite_index_single, composite_index_array, all_composite_indices
        generic :: composite_index => composite_index_single, composite_index_array, all_composite_indices

      procedure :: multi_index_single, multi_index_array, all_multi_indices
        generic :: multi_index => multi_index_single, multi_index_array, all_multi_indices
    end type regular_mesh_type

    ! Overload default constructor
    !> Setup an instance of `[[regular_mesh_type]]`.
    interface regular_mesh_type
      procedure :: setup_regular_mesh
    end interface regular_mesh_type

    contains

    !> Setup an instance of `[[regular_mesh_type]]`.
    function setup_regular_mesh(sampling, multi_index_first) result(this)
      !> Number of grid points per dimension.
      integer, intent(in) :: sampling(3)
      !> Multi index, assigned to `composite_index = 1` for indentity mapping.
      integer, intent(in) :: multi_index_first(3)

      type(regular_mesh_type) :: this

      CALL_ASSERT(all(sampling > 0), 'sampling <= 0.')
      CALL_ASSERT(all(multi_index_first > 0), 'multi_index_first <= 0.')
      CALL_ASSERT(all(multi_index_first <= sampling), 'multi_index_first > sampling.')

      this%sampling = sampling
      this%multi_index_first = multi_index_first
    end function setup_regular_mesh

    !> Return the number of points in the mesh.
    integer function number_of_points(this)
      class(regular_mesh_type), intent(in) :: this

      number_of_points = product(this%sampling)
    end function number_of_points

    !---------------------------------------------------------------------------------------------
    ! Composite index access
    !---------------------------------------------------------------------------------------------

    !> Map `multi_index` to `composite_index`.
    function composite_index_single(this, multi_index) result(composite_index)
      !> Mesh
      class(regular_mesh_type), intent(in) :: this
      !> Integer coordinate
      integer, intent(in) :: multi_index(3)

      integer :: composite_index

      integer :: multi_index_local(3)

      multi_index_local = mod1(multi_index - this%multi_index_first + 1, this%sampling)
      composite_index = indices_to_composite_index(multi_index_local, this%sampling)
    end function composite_index_single

    !> Map an array of `multi_indices` to array of `composite_indices`.
    function composite_index_array(this, multi_indices) result(composite_indices)
      class(regular_mesh_type), intent(in) :: this
      integer, intent(in) :: multi_indices(:, :)

      integer, allocatable :: composite_indices(:)

      integer :: number_of_points, i

      CALL_ASSERT(size(multi_indices, 1) == 3, 'size(multi_indices, 1) /= 3.')

      number_of_points = size(multi_indices, 2)
      composite_indices = [(this%composite_index(multi_indices(:, i)), i = 1, number_of_points)]
    end function composite_index_array

    !> Return all composite indices of the grid.
    function all_composite_indices(this) result(composite_indices)
      class(regular_mesh_type), intent(in) :: this

      integer, allocatable :: composite_indices(:)

      composite_indices = mesh_1d(1, this%number_of_points())
    end function all_composite_indices

    !---------------------------------------------------------------------------------------------
    ! Multi index access
    !---------------------------------------------------------------------------------------------

    !> Map `composite_index` to `multi_index`.
    function multi_index_single(this, composite_index) result(multi_index)
      class(regular_mesh_type), intent(in) :: this
      integer, intent(in) :: composite_index

      integer :: multi_index(3)

      integer :: composite_index_local

      composite_index_local = mod1(composite_index, this%number_of_points())
      call composite_index_to_indices(composite_index_local, this%sampling, multi_index)
      multi_index = mod1(multi_index + this%multi_index_first - 1, this%sampling)
    end function multi_index_single

    !> Map array of `composite_indices` to array of `multi_indices`.
    function multi_index_array(this, composite_indices) result(multi_indices)
      class(regular_mesh_type), intent(in) :: this
      integer, intent(in) :: composite_indices(:)

      integer, allocatable :: multi_indices(:, :)

      integer :: number_of_points, i

      number_of_points = size(composite_indices)
      allocate(multi_indices(3, number_of_points))

      do i=1, number_of_points
        multi_indices(:, i) = this%multi_index(composite_indices(i))
      end do
    end function multi_index_array

    !> Return all multi indices of the grid
    function all_multi_indices(this) result(multi_indices)
      class(regular_mesh_type), intent(in) :: this

      integer, allocatable :: multi_indices(:, :)

      multi_indices = this%multi_index(mesh_1d(1, this%number_of_points()))
    end function all_multi_indices

end module regular_mesh