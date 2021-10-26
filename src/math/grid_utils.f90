!> Tools that are useful for defining grids or working with grids.
module grid_utils
  use iso_fortran_env, only: error_unit

  use precision, only: sp, dp
  use constants, only: pi, zi
  use asserts, only: assert
  use math_utils, only: kronecker_product, mod1
  
  implicit none
  
  private
  public :: linspace, &
            grid_3d, &
            phase, &
            fft_frequencies
            

  !> Generate a vector of evenly spaced number over a given interval. 
  !> The spacing can be either defined directly or by a number of points.
  interface linspace
    module procedure :: linspace_spacing, linspace_number_of_points
  end interface linspace

  interface grid_3d 
    module procedure :: cubic_grid_3d, regular_cubic_grid_3d
  end interface grid_3d

  interface phase
    module procedure :: phase_single_point, phase_array
  end interface phase
  
contains

  !> Return a vector of evenly-spaced numbers over a given interval.
  !>
  !> The start point is always included in the grid, whereas the end point 
  !> is never included. If start == end, a single point is returned, 
  !> irrespective of spacing.
  !>
  !> If start > end, the grid is returned in descending from start to end.
  !> For example:
  !>
  !>  [5, 4, 3, 2, 1] = linspace_spacing(5, 0, 1) 
  !>
  function linspace_spacing(start, end, spacing) result(grid) 
    !> Start of the grid
    real(dp), intent(in) :: start
    !> End of the grid: This number is never contained in the grid.
    real(dp), intent(in) ::  end
    !> Spacing between the points
    real(dp), intent(in) :: spacing
    !> Grid
    real(dp), allocatable :: grid(:)

    integer :: N, i, sign

    if (start == end) then
      allocate(grid(1), source = start)
      return 
    end if 

    call assert(spacing > 0._dp, & 
         message='linspace: Spacing must be larger than zero.')

    call assert(abs(start - end) >= spacing, &
         message = "linspace: Requires spacing <= |start - end|")

    N = floor(abs(end - start) / spacing) + 1
    allocate(grid(N)) 
 
    sign = 1 
    if (end < start) then
      sign = -1
    endif
 
    do i = 1, N
      grid(i) = start + sign * (i - 1) * abs(spacing)
    end do

  end function linspace_spacing


  !> Generate a linearly-spaced grid between 'start' and 'end'. 
  !> 
  !> The number of grid points, N, must be > 0.
  !> Default is to include the endpoint
  function linspace_number_of_points(start, end, N, endpoint) result(grid)
    !> Start of the grid
    real(dp), intent(in) :: start
    !> End of the grid
    real(dp), intent(in) :: end 
    !> Number of points
    integer, intent(in) :: N
    !> If true, include the endpoint 'end' in the grid 
    logical, intent(in), optional :: endpoint
    !> Linear grid
    real(dp), allocatable :: grid (:), grid_(:)

    logical :: include_endpoint
    real(dp) :: spacing
    
    call assert(N > 0, &
                  message = 'linspace: The number of grid points must not be zero or smaller (N>0).')

    include_endpoint = .true.
    if (present(endpoint)) include_endpoint = endpoint

    allocate(grid(N))
    
    if (N == 1) then
      grid = [start]

    else
      if (include_endpoint) then
        spacing = abs(end - start) / real(N - 1, kind=dp)
        grid(1: N-1) = linspace_spacing(start, end, spacing)  
        grid(N) = end

      else   
        spacing = abs(end - start) / real(N, kind=dp)
        grid_ = linspace_spacing(start, end, spacing)
        grid = grid_(:N)

      end if
    end if

  end function linspace_number_of_points


  !> Generate a cubic, 3D grid.
  !>  
  !> Sampling be specified per dimension by N.
  !>
  !> The grid is defined by the number of grid points per dimension and a cuboid.
  !> The cuboid is defined by the points 'start' and 'end', which span the diagonal 
  !> of the cuboid.
  !> 
  !> The defaults for the start and of the cubic are (\(0, 0, 0)\) and (\(1, 1, 1)\),
  !> respectively.
  !>
  !> The grid containes by default the end point. If this behavior is not wished, 
  !> the end point can be left out by setting endpoint to .false.
  function cubic_grid_3d(N, start, end, endpoint) result(grid)
    !> Number of grid points per dimension
    integer, intent(in) :: N(3)
    !> Define the start point of the grid in each dimension
    real(dp), intent(in), optional :: start(3)
    !> Define the end point of the grid in each dimension
    real(dp), intent(in), optional :: end(3)
    !> If true, include the endpoint 'end' in the grid 
    logical, intent(in), optional :: endpoint
    !> Cubic, 3D grid 
    real(dp), allocatable :: grid(:,:)

    !> Default origin for the grid
    real(dp), dimension(3), parameter :: default_start = [0._dp, 0._dp, 0._dp]
    !> Default end point for the grid 
    real(dp), dimension(3), parameter :: default_end   = [1._dp, 1._dp, 1._dp]

    real(dp) :: start_(3), end_(3), ones_x(N(1)), ones_y(N(2)), ones_z(N(3)) 
    real(dp), allocatable :: grid_1D_x(:), grid_1D_y(:), grid_1D_z(:)
      
    logical :: include_endpoint
    
    include_endpoint = .true.
    if (present(endpoint)) include_endpoint = endpoint

    ! allocate grid
    allocate(grid(3, product(N)))

    ! set start and end for constructing the cuboid
    start_ = default_start
    end_ = default_end 
    if (present(start)) start_ = start
    if (present(end)) end_ = end

    ! Generate 1d grids for each direction
    grid_1D_x = linspace(start_(1), end_(1), N(1), endpoint = include_endpoint)
    grid_1D_y = linspace(start_(2), end_(2), N(2), endpoint = include_endpoint)
    grid_1D_z = linspace(start_(3), end_(3), N(3), endpoint = include_endpoint)
            
    ! put 1d grids together
    ones_x = 1.0_dp
    ones_y = 1.0_dp
    ones_z = 1.0_dp

    grid(1,:) = kronecker_product( ones_z,    ones_y,    grid_1D_x)
    grid(2,:) = kronecker_product( ones_z,    grid_1D_y, ones_x)
    grid(3,:) = kronecker_product( grid_1D_z, ones_y,    ones_x)        
  end function cubic_grid_3d


  !> Generate a regular cubic grid.
  !>  
  !> See [[cubic_grid_3d]]
  function regular_cubic_grid_3d(N, start, end, endpoint) result(grid)
    !> Number of grid points per dimension
    integer, intent(in) :: N
    !> Define the start point of the grid in each dimension
    real(dp), intent(in), optional :: start(3)
    !> Define the end point of the grid in each dimension
    real(dp), intent(in), optional :: end(3)
    !> If true, include the endpoint 'end' in the grid 
    logical, intent(in), optional :: endpoint
    !> Regular cubic grid 
    real(dp), allocatable :: grid(:, :)
    logical :: include_endpoint
    
    include_endpoint = .true.
    if (present(endpoint)) include_endpoint = endpoint

    grid = cubic_grid_3d([N, N, N], start, end)
  end function

  !> Generate the Bloch phase for a point in the unit cell and a k-point:
  !> \[ 
  !>     \text{phase}(\mathbf r, \mathbf k) = 
  !>        e^{-2 \pi i \mathbf k \cdot \mathbf r}, 
  !> \]
  !> where
  !> \[ 
  !>       \mathbf r_i \in \text{unit cell},  \mathbf k \in \text{BZ}. 
  !> \]
  complex(dp) function phase_single_point(r, k)
    !> r-point
    real(dp), intent(in) :: r(3)
    !> k-point
    real(dp), intent(in) :: k(3)

    phase_single_point = exp(-2 * zi * pi * sum(k * r))
  end function phase_single_point

  !> Generate the complex phase on a grid for a k-point and store it to an array in the correct order:
  !> \[ 
  !>     \text{phase}(\mathbf r_i, \mathbf k) = 
  !>        e^{-2 \pi i \mathbf k \cdot \mathbf r_i}, 
  !> \]
  !> where
  !> \[ 
  !>       \mathbf r_i \in \text{grid},  \mathbf k \in \text{BZ}. 
  !> \]
  function phase_array(r_array, k) result(phase)
    !> Array of r-points
    real(dp), intent(in) :: r_array(:,:)
    !> k_point
    real(dp), intent(in) :: k(3)
    !> Phase
    complex(dp), allocatable :: phase(:)

    integer :: i

    call assert(size(r_array, 1) == 3, &
                'phase_array: First dimension of r_array needs to be 3.')
    
    allocate(phase(size(r_array, 2)))

    do i=1, size(r_array, 2)
      phase(i) = phase_single_point(r_array(:, i), k)
    end do
  end function phase_array

  !> Calculate Fourier frequencies of a function on a three dimensional grid 
  !> in the order as they appear in the 3D FFT for the given grid.
  function fft_frequencies(N) result(frequencies)
    implicit none
    !> Number of grid points per dimension
    integer, intent(in) :: N(3)
    !> Fourier frequncies
    real(dp) :: frequencies(3, product(N))
    ! local variables
    integer :: pos_max(3), neg_min(3)
    real(dp) :: ones1(N(1)), ones2(N(2)), ones3(N(3)), &
                f1(N(1)), f2(N(2)), f3(N(3))
    integer :: i
      
    do i = 1, 3
        pos_max(i) =  int(floor(real(N(i), kind=sp) / 2.0_sp)) - (1 - modulo(N(i), 2))
        neg_min(i) = -int(floor(real(N(i), kind=sp) / 2.0_sp))
    end do  

    ones1 = 1.0_dp
    f1(1:pos_max(1)+1) = [(i, i=0, pos_max(1))]
    f1(pos_max(1)+2:N(1)) = [(i, i=neg_min(1), -1)]

    ones2 = 1.0_dp
    f2(1:pos_max(2)+1) = [(i, i=0, pos_max(2))]
    f2(pos_max(2)+2:N(2)) = [(i, i=neg_min(2), -1)]

    ones3 = 1.0_dp
    f3(1:pos_max(3)+1) = [(i, i=0, pos_max(3))]
    f3(pos_max(3)+2:N(3)) = [(i, i=neg_min(3), -1)]

    frequencies(1, :) = kronecker_product(ones3, ones2, f1)
    frequencies(2, :) = kronecker_product(ones3, f2,    ones1)
    frequencies(3, :) = kronecker_product(f3,    ones2, ones1)
  end function fft_frequencies

end module grid_utils
