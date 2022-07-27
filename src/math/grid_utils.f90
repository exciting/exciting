!> Tools that are useful for defining grids or working with grids.
module grid_utils
  use precision, only: sp, dp
  use constants, only: pi, zi
  use asserts, only: assert
  use math_utils, only: kronecker_product, mod1
  use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices

  implicit none
  
  private
  public :: mesh_1d, &
            linspace, &
            concatenate, &
            grid_3d, &
            phase, &
            fft_frequencies, &
            n_grid_diff, &
            partial_grid, &
            flattend_map

  real(dp), parameter :: default_offset(3) = [0._dp, 0._dp, 0._dp]

  !> Generate a vector of evenly spaced number over a given interval.
  !> The spacing can be either defined directly or by a number of coords.
  interface linspace
    module procedure :: linspace_spacing, linspace_number_of_points, linspace_interval
  end interface linspace

  !> Concatanate two 1-rank arrays:
  !>
  !> `concatenate([v1, ..., v_n], [w1, ..., w_m]) -> [v1, ..., v_n, w_1, ..., w_m]`
  interface concatenate
    module procedure :: concatenate_integer, concatenate_real_dp, concatenate_complex_dp
  end interface concatenate

  !> Generate a regular grid with the last index changing the fastest
  interface grid_3d
    module procedure :: cubic_grid_3d, regular_cubic_grid_3d
  end interface grid_3d

  !> Calculate the complex phase for a given k-point and a r-point or
  !> an array of r-points
  interface phase
    module procedure :: phase_single_point, phase_array
  end interface phase
  
contains

  !> Generate an array of ascending or descending evenly distributed integers.
  !> The first element is alway `start` and the last element `mesh_1d_last <= end` (for `end < start`:
  !> `mesh_1d_last <= end`).
  !>
  !> `mesh_1d_last == end` when `(start - end) / spacing` is an integer number.
  function mesh_1d(start, end, spacing) result(range)
    !> First element of the range. If end is not given, this is interpreted as length of the mesh 
    !> and the first element is set to 1
    integer, intent(in) :: start
    !> Last element of the range. Must be different from start. If `end < start`, a descending array is returned.
    integer, intent(in) :: end
    !> Spacing between the elements. The routine expects `spacing > 0`. Default is 1.
    integer, intent(in), optional :: spacing

    integer, allocatable :: range(:)

    integer :: i, spacing_local, dx, N
    external :: sign

    spacing_local = 1
    if (present(spacing)) spacing_local = spacing

    call assert(start /= end, 'start == end.')
    call assert(spacing_local > 0, 'spacing <= 0.')

    dx = end - start
    spacing_local = isign(spacing_local, dx)
    N = int(dx / spacing_local) + 1

    allocate(range(N))
    do i = 1, N
       range(i) = start + (i - 1) * spacing_local
    end do
  end function mesh_1d


  !> Generate an array of lenght N_tot of ascending distributed intergers
  !> from 1 to N_tot.
  !> This is generating the identiy map.
  function mesh_1d_N_points(N_points) result(range_out)
    !> First element of the range
    integer, intent(in) :: N_points

    integer, allocatable :: range_out(:)

    integer :: i

    call assert(N_points > 0, 'N_points <= 0.')
    allocate(range_out(N_points))

    do i=1, N_points
      range_out(i) = i
    end do
  end function mesh_1d_N_points

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
  function linspace_number_of_points(start, end, N, include_end) result(grid)
    !> Start of the grid
    real(dp), intent(in) :: start
    !> End of the grid
    real(dp), intent(in) :: end 
    !> Number of points
    integer, intent(in) :: N
    !> If true, `end` is the last element of the grid and `spacing = (end - start) / (N - 1)`, else
    !> the last element is given by `end - spacing`, where`spacing = (end - start) / N`.
    !> 
    !> Default is true.
    logical, intent(in), optional :: include_end
    !> Linear grid
    real(dp), allocatable :: grid (:), grid_(:)

    logical :: include_end_local
    real(dp) :: spacing
    
    call assert(N > 0, 'The number of grid points (N) is <= 0.')

    include_end_local = .true.
    if (present(include_end)) include_end_local = include_end

    allocate(grid(N))
    
    if (N == 1) then
      grid = [start]

    else
      if (include_end_local) then
        spacing = abs(end - start) / real(N - 1, kind=dp)
        grid_ = linspace_spacing(start, end, spacing)
        grid(1: N-1) = grid_(1: N-1)
        grid(N) = end

      else   
        spacing = abs(end - start) / real(N, kind=dp)
        grid_ = linspace_spacing(start, end, spacing)
        grid = grid_(:N)

      end if
    end if
  end function linspace_number_of_points

  !> Generate a linearly-spaced grid in the interval `interval = [start, end]`.
  function linspace_interval(interval, N, include_end) result(grid)
    !> Interval to define the grid
    real(dp), intent(in) :: interval(2)
    !> Number of coords
    integer, intent(in) :: N
    !> Influences the spacing between the points. If the include_end is included, `spacing = (end - start) / (N - 1)`,
    !> else `spacing = (end - start) / N`.
    logical, intent(in), optional :: include_end
    !> Linear grid
    real(dp), allocatable :: grid (:)
    !> If true, `end` is the last element of the grid and `spacing = (end - start) / (N - 1)`, else
    !> the last element is given by `end - spacing`, where`spacing = (end - start) / N`.
    !> 
    !> Default is true.
    logical :: include_end_local

    include_end_local = .true.
    if (present(include_end)) include_end_local = include_end

    grid = linspace_number_of_points(interval(1), interval(2), N, include_end_local)
  end function linspace_interval

  !> Concatenate two integer 1-rank arrays.
  function concatenate_integer(v1, v2) result(v_out)
    !> Vectors to be concatenated.
    integer, intent(in) :: v1(:), v2(:)

    integer, allocatable :: v_out(:)

    integer :: len_1, len_2

    len_1 = size(v1)
    len_2 = size(v2)

    call assert(len_1 > 0, 'size(v1) == 0.')
    call assert(len_2 > 0, 'size(v2) == 0.')

    allocate(v_out(len_1 + len_2))

    v_out(: len_1) = v1
    v_out(len_1 + 1 :) = v2
  end function concatenate_integer

  !> Concatenate two double real 1-rank arrays.
  function concatenate_real_dp(v1, v2) result(v_out)
    !> Vectors to be concatenated.
    real(dp), intent(in) :: v1(:), v2(:)

    real(dp), allocatable :: v_out(:)

    integer :: len_1, len_2

    len_1 = size(v1)
    len_2 = size(v2)

    call assert(len_1 > 0, 'size(v1) == 0.')
    call assert(len_2 > 0, 'size(v2) == 0.')

    allocate(v_out(len_1 + len_2))

    v_out(: len_1) = v1
    v_out(len_1 + 1 :) = v2
  end function concatenate_real_dp

  !> Concatenate two double complex 1-rank arrays.
  function concatenate_complex_dp(v1, v2) result(v_out)
    !> Vectors to be concatenated.
    complex(dp), intent(in) :: v1(:), v2(:)

    complex(dp), allocatable :: v_out(:)

    integer :: len_1, len_2

    len_1 = size(v1)
    len_2 = size(v2)

    call assert(len_1 > 0, 'size(v1) == 0.')
    call assert(len_2 > 0, 'size(v2) == 0.')

    allocate(v_out(len_1 + len_2))

    v_out(: len_1) = v1
    v_out(len_1 + 1 :) = v2

    v_out = reshape([v1, v2], [len_1 + len_2])
  end function concatenate_complex_dp


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
  !> the end point can be left out by setting include_end to .false.
  function cubic_grid_3d(N, start, end, include_end) result(grid)
    !> Number of grid points per dimension
    integer, intent(in) :: N(3)
    !> Define the start point of the grid in each dimension
    real(dp), intent(in), optional :: start(3)
    !> Define the end point of the grid in each dimension
    real(dp), intent(in), optional :: end(3)
    !> If true, include 'end' in the grid (see [[linspace_number_of_points]]).
    logical, intent(in), optional :: include_end
    !> Cubic, 3D grid 
    real(dp), allocatable :: grid(:,:)

    !> Default origin for the grid
    real(dp), dimension(3), parameter :: default_start = [0._dp, 0._dp, 0._dp]
    !> Default end point for the grid
    real(dp), dimension(3), parameter :: default_end   = [1._dp, 1._dp, 1._dp]

    real(dp) :: start_(3), end_(3), ones_x(N(1)), ones_y(N(2)), ones_z(N(3)) 
    real(dp), allocatable :: grid_1D_x(:), grid_1D_y(:), grid_1D_z(:)
      
    logical :: include_end_local
    
    include_end_local = .true.
    if (present(include_end)) include_end_local = include_end

    ! allocate grid
    allocate(grid(3, product(N)))

    ! set start and end for constructing the cuboid
    start_ = default_start
    end_ = default_end 
    if (present(start)) start_ = start
    if (present(end)) end_ = end

    ! Generate 1d grids for each direction
    grid_1D_x = linspace(start_(1), end_(1), N(1), include_end = include_end_local)
    grid_1D_y = linspace(start_(2), end_(2), N(2), include_end = include_end_local)
    grid_1D_z = linspace(start_(3), end_(3), N(3), include_end = include_end_local)
            
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
  function regular_cubic_grid_3d(N, start, end, include_end) result(grid)
    !> Number of grid points per dimension
    integer, intent(in) :: N
    !> Define the start point of the grid in each dimension
    real(dp), intent(in), optional :: start(3)
    !> Define the end point of the grid in each dimension
    real(dp), intent(in), optional :: end(3)
    !> If true, include 'end' in the grid (see [[linspace_number_of_points]]).
    logical, intent(in), optional :: include_end
    !> Regular cubic grid 
    real(dp), allocatable :: grid(:, :)
    logical :: include_end_local
    
    include_end_local = .true.
    if (present(include_end)) include_end_local = include_end

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

    call assert(size(r_array, 1) == 3, 'First dimension of r_array needs to be 3.')
    
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


  !> Calculate the number of grid points per dimension of the difference grid with the given the number of grid
  !> points per dimension of a regular origin grid. In each dimension with more then one grid point, the number of grid
  !> points is doubled because the difference grid take to into account \( \mathbf{p}_1 - \mathf{p}_2\) and
  !> \( \mathbf{p}_2 - \mathf{p}_1\) where \(\mathbf{p}_1, \mathbf{p}_2\) are points of the origin grid. Example:
  !> \[
  !>   \begin{split}
  !>     2 3 4 \rightarrow 4 6 8
  !>     1 2 3 \rightarrow 1 4 6
  !>   \end{split}
  !> \]
  function n_grid_diff(N_grid_in) result(N_grid_out)
    !> Number of grid points per dimension
    integer, intent(in) :: N_grid_in(3)

    integer :: N_grid_out(3)

    integer :: i

    call assert(all(N_grid_in > 0),  'Number of k-points per dimension is zero.')
    N_grid_out = 2 * N_grid_in
    do i = 1, 3
      if (N_grid_in(i) == 1) N_grid_out(i) = 1
    end do
  end function n_grid_diff

  !> Calculate for a grid that is included in a denser grid, the indices of the grid points
  !> in the denser grid. The number of grid points of the grid must divide the number of
  !> grid points in the denser grid to be included.
  function partial_grid(N_grid_dense, N_grid) result(map)
    !> Number of grid points per dimension of the bigger grid
    integer, intent(in) :: N_grid_dense(3)
    !> Number of grid points per dimension of the smaller grid
    integer, intent(in) :: N_grid(3)

    integer, allocatable :: map(:)

    integer :: divider(3), i, index3d(3)

    call assert(all(mod(N_grid_dense, N_grid) /= 0) , 'N_grid_dense modulus N_grid is not zero for all dimensions.')

    divider = N_grid_dense / N_grid
    allocate(map(product(N_grid)))

    do i = 1, product(N_grid)
      call composite_index_to_indices(i, N_grid, index3d)
      map(i) = indices_to_composite_index(index3d * divider, N_grid_dense)
    end do
  end function partial_grid

  !> Given the shapes of two `M`-rank arrays `shape_A` and `shape_B` with `all(shape_B <= shape_A)` is true and
  !> `A(1 : shape_B(1), ..., 1 : shape_B(M)) = B, find the index map `map_out` of length `product(shape_B)` such that 
  !> 'A_flat(map_out) = B_flat`. Where `A_flat = reshape(A, [product(shape_A)])` and `B_flat = reshape(B, [product(shape_B)])`.
  function flattend_map(shape_A, shape_B) result(map_out)
    !> Shape of `A`
    integer, intent(in) :: shape_A(:)
    !> Shape of `B`
    integer, intent(in) :: shape_B(:)

    integer, allocatable :: map_out(:)

    integer :: i, rank, size_B
    integer, allocatable :: multi_index(:)


    call assert(size(shape_B) == size(shape_A), 'shape_B and shape_A have not the same size.')
    call assert(all(shape_A >= shape_B), 'shape_A < shape_B.')

    size_B = product(shape_B)
    rank = size(shape_B)

    allocate(multi_index(rank))
    allocate(map_out(size_B))

    do i = 1, size_B
      ! Calculate the multi index of B from the index for B_flat
      call composite_index_to_indices(i, rank, shape_B, multi_index)
      ! Calculate the index of A_flat from the multi index of B
      map_out(i) = indices_to_composite_index(multi_index, shape_A)
    end do
  end function flattend_map

end module grid_utils
