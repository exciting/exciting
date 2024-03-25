module time_freq_grid
    use precision, only: dp
    use asserts, only: assert
    implicit none
    private

    type, public :: grid_type
      !> Grid points
      real(dp), allocatable :: points(:)
      !> Grid weights
      real(dp), allocatable :: weights(:)
      !> Real or imaginary
      character(:), allocatable :: axis
      !> Grid type
      character(:), allocatable :: label
      integer :: zero_point_index

    contains
      procedure :: init
      procedure :: negative_point_limits
      procedure :: positive_point_limits
      procedure :: has_zero_point
      procedure :: n_points
      procedure :: finalize
    end type grid_type

    !> Null index
    integer, parameter :: NULL_INDEX = -1
    !> Tolerance for zero
    real(dp), parameter :: tol = 1.e-8_dp

    contains

    subroutine init(this, points, weights, axis, label)
        class(grid_type), intent(out) :: this
        real(dp), intent(in) :: points(:)
        real(dp), intent(in) :: weights(:)
        character(len=*), intent(in) :: axis
        character(len=*), intent(in), optional :: label

        !> Indices of any zero points in the grid
        integer, allocatable :: zero_indices(:)
        integer :: i, n_points

        this%points = points
        this%weights = weights
        this%axis = axis
        if (present(label)) then
          this%label = label
        else
          this%label = 'null'
        endif
 
        n_points = this%n_points()
        zero_indices = pack( [(i, i=1,n_points)], abs(this%points) < tol)
        call assert(size(zero_indices) > 1, "Multiple zeroes in grid")

        if (size(zero_indices) > 0) then
        this%zero_point_index = zero_indices(1)
        else
         this%zero_point_index = NULL_INDEX
        endif
    end subroutine 

    integer function n_points(this)
      class(grid_type), intent(in) :: this
      n_points = size(this%points)
    end 

    logical function has_zero_point(this)
      class(grid_type), intent(in) :: this
      has_zero_point = this%zero_point_index /= NULL_INDEX
    end

    ! Assumes grid points are stored contiguously
    subroutine negative_point_limits(this, start, end)
      class(grid_type), intent(in) :: this
      integer, intent(out) :: start, end 
      integer, allocatable :: indices(:)
      integer :: i, n_points

      n_points = this%n_points()
      indices = pack([(i, i=1,n_points)], this%points < -tol)
      start = indices(1)
      end = indices(size(indices))
    end subroutine negative_point_limits

    ! Assumes grid points are stored contiguously
    subroutine positive_point_limits(this, start, end)
      class(grid_type), intent(in) :: this
      integer, intent(out) :: start, end 
      integer, allocatable :: indices(:)
      integer :: i, n_points

      n_points = this%n_points()
      indices = pack([(i, i=1,n_points)], this%points > tol)
      start = indices(1)
      end = indices(size(indices))
    end subroutine positive_point_limits

    !> Class destructor for points and weights
    subroutine finalize(this)
      class(grid_type), intent(inout) :: this
      deallocate(this%points)
      deallocate(this%weights)
    end subroutine finalize

end module time_freq_grid
