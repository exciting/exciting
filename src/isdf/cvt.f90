module cvt_utils
  use precision, only: dp
  use asserts, only: assert
  use math_utils, only: calculate_all_vector_distances, all_close
  use grid_utils, only: mesh_1d

  implicit none

  private
  public :: cvt

  interface cvt
    procedure :: cvt_interpolation_point_indices, cvt_procedure
  end interface cvt

  contains

  
  subroutine cvt_interpolation_point_indices(points, weights, convergence_tolerance, steps, interpolation_point_indices, has_converged, deviation)
    !> Coordinates of the points that are assigned to the clusters.
    real(dp), contiguous, intent(in) :: points(:, :) ! 3 x n_points
    !> Weights corresponding to the points.
    real(dp), contiguous, intent(in) :: weights(:) ! n_points
    !> Tolerance for deciding if the cluster centroids are element wise the same as in the previous step.
    real(dp), intent(in) :: convergence_tolerance
    !> On input, maximum steps to perform. If reached, the procedure returns the last calculated cluster centroids.
    !> On output, number of steps, acutally performed.
    integer, intent(inout) :: steps
    !> On input, indices to elements of `points` that are used as initial cluster centroids.
    !> On output, indices to elements of `points` that are closest to the cluster coordinates as calculated by `[[cvt_procedure]]`.
    integer, intent(inout) :: interpolation_point_indices(:)
    !> Flag that determines if `[[cluster_centroids]]` is converged or not after return.
    !> It is converged if `[[cluster_centroids]]` calculated in the current step are element wise the same as
    !> as in the previous step, within the `[[convergence_tolerance]]`, before reaching `[[maximum_steps]]`. 
    logical, intent(out) :: has_converged
    !> Mean deviation of the cluster centroids as calculated by `[[cvt_procedure]]` and `[[interpolation_point_indices]]` 
    !> point to on output.
    real(dp), intent(out) :: deviation

    integer :: n_points, n_interpolation_points
    real(dp), allocatable :: cluster_centroids(:, :), distance_matrix(:, :)

    call assert(size(points, 1) == 3, &
            'Columns of points do not contain 3-d cooridnates: size(points, 1) \= 3.')

    call assert(size(points, 2) >= size(interpolation_point_indices), &
            'The number of requested interpolation points is larger then the number of given points: &
                    size(points, 2) < size(interpolation_point_indices)')

    call assert(all(interpolation_point_indices <= size(points, 2)), &
            'Some elements in interpolation_point_indices are larger then the total number of points: &
                    any(interpolation_point_indices > size(points, 2))')

    call assert(all(interpolation_point_indices > 0), &
            'Some elements in interpolation_point_indices are less smaller then zero: &
                    any(interpolation_point_indices <= 0)')

    ! TODO: assert that interpolation_point_indices is an injective map.

    n_points = size(points, 2)                
    n_interpolation_points = size(interpolation_point_indices)
    allocate(distance_matrix(n_points, n_interpolation_points))

    ! Run CVT procedure
    cluster_centroids = points(:, interpolation_point_indices)
    call cvt(points, weights, convergence_tolerance, steps, cluster_centroids, has_converged)

    ! Map cluster centroids to the closest points
    call calculate_all_vector_distances(points, cluster_centroids, distance_matrix)
    interpolation_point_indices = minloc(distance_matrix, dim = 1)

    ! Calculate deviation
    deviation = sum(norm2(points(:, interpolation_point_indices) - cluster_centroids, dim=1)) / n_interpolation_points

    deallocate(cluster_centroids, distance_matrix)
  end subroutine

  subroutine cvt_procedure(points, weights, convergence_tolerance, steps, cluster_centroids, has_converged)
    !> Coordinates of the points that are assigned to the clusters.
    real(dp), contiguous, intent(in) :: points(:, :) ! 3 x n_points
    !> Weights corresponding to the points.
    real(dp), contiguous, intent(in) :: weights(:) ! n_points
    !> Tolerance for deciding if the cluster centroids are element wise the same as in the previous step.
    real(dp), intent(in) :: convergence_tolerance
    !> On input, maximum steps to perform. If reached, the procedure returns the last calculated cluster centroids.
    !> On output, number of steps, acutally performed.
    integer, intent(inout) :: steps
    !> On input, cartesian cooridinate of the centroids of the inital clusters, 
    !> on output cartesian coordinates of the centroids of the optimized clusters.
    real(dp), contiguous, intent(inout) :: cluster_centroids(:, :) ! 3 x n_clusters
    !> Flag that determines if `[[cluster_centroids]]` is converged or not after return.
    !> It is converged if `[[cluster_centroids]]` calculated in the current step are element wise the same as
    !> as in the previous step, within the `[[convergence_tolerance]]`, before reaching `[[steps]]`. 
    logical, intent(out) :: has_converged
    

    integer :: step, n_clusters, n_points, maximum_steps
    integer, allocatable :: map_to_clusters(:)
    real(dp), allocatable :: old_cluster_centroids(:, :), weighted_points(:, :), distance_matrix(:, :)
    
    call assert(size(points, 1) == 3, &
            'Columns of points are not 3-dim points: size(points, 1) /= 3.')

    call assert(size(points, 2) == size(weights), &
            'Not the same number of points and weights: size(points, 2) /= size(weights).')

    call assert(size(cluster_centroids, 1) == 3, &
            'Columns of cluster_centroids are not 3-dim points: size(cluster_centroids, 1) /= 3.')

    call assert(size(points, 2) >= size(cluster_centroids, 2), &
            'Number of points is smaller then number of clusters: size(points, 2) < size(cluster_centroids, 2).')

    call assert(steps > 0, 'steps <= 0 on input.')

    n_points = size(points, 2)
    n_clusters = size(cluster_centroids, 2)

    allocate(distance_matrix(n_clusters, n_points))
    allocate(map_to_clusters(n_points))

    weighted_points = points * spread(weights, dim=1, ncopies=3)
    maximum_steps = steps
    step = 0
    has_converged = .false.
    do while(.not. has_converged)

      call calculate_all_vector_distances(cluster_centroids, points, distance_matrix)
      call assign_points_to_cluster(distance_matrix, map_to_clusters)

      ! Assign via pointer to avoid copying?
      old_cluster_centroids = cluster_centroids
      call update_cluster_cenctroids(weighted_points, weights, map_to_clusters, cluster_centroids)

      has_converged = all_close(cluster_centroids, old_cluster_centroids, tol=convergence_tolerance)

      if(step == maximum_steps) exit
      step = step + 1
    end do

    steps = step

    deallocate(map_to_clusters, old_cluster_centroids, weighted_points)
  end subroutine cvt_procedure

  subroutine update_cluster_cenctroids(weigthed_points, weights, map_to_clusters, cluster_centroids)
    !> Column-wise cartesian point cooridinates, multiplied by their weights.
    real(dp), contiguous, intent(in) :: weigthed_points(:, :) ! 3 x n_points
    !> Weights corresponding to the points.
    real(dp), contiguous, intent(in) :: weights(:) ! n_points
    !> Index map from point to cluster. `map_to_cluster(5) = 3` assigns `point(:, 5)` to 
    !> the cluster corresponding to `cluster_centroids(:, 3)`.
    integer, contiguous, intent(in) :: map_to_clusters(:) ! n_points
    !> Cartesian cluster centroid coordinates.
    real(dp), contiguous, intent(out) :: cluster_centroids(:, :) ! 3 x n_clusters
    
    integer :: i_cluster, n_clusters, n_points, i_point, ind

    integer, allocatable :: indices(:), indices_in_cluster(:)

    real(dp) :: cluster_weight, cluster_centroid(3)

    n_clusters = size(cluster_centroids, 2)
    n_points = size(weights)

    indices = mesh_1d(1, n_points)

    do  i_cluster = 1, n_clusters

      indices_in_cluster = pack(indices, map_to_clusters == i_cluster)

      cluster_centroid = 0._dp
      cluster_weight = 0._dp

      ! Start the parallel region in the outer loop to reduce fork:merge?
      !$OMP parallel do reduction(+: cluster_centroid, cluster_weight) shared(indices_in_cluster, weigthed_points) private(i_point)
      do ind=1, size(indices_in_cluster)
        i_point = indices_in_cluster(ind)

        cluster_centroid = cluster_centroid + weigthed_points(:, i_point)
        cluster_weight = cluster_weight + weights(i_point)
      end do
      !$OMP end parallel do

      cluster_centroids(:, i_cluster) = cluster_centroid / cluster_weight

    end do

  end subroutine update_cluster_cenctroids

  subroutine assign_points_to_cluster(distance_matrix, map_to_cluster)
    !> Distance matrix. Each row contains the distances of a point to all cluster centroids.
    real(dp), contiguous, intent(in) :: distance_matrix(:, :) ! n_clusters x n_points
    !> Index map from point to cluster. `map_to_cluster(5) = 3` assigns `point(:, 5)` to 
    !> the cluster corresponding to `cluster_centroids(:, 3)`.
    integer, contiguous, intent(out) :: map_to_cluster(:) ! n_points

    integer :: i_point
  
    ! single for opening parallel outside 
    if(size(distance_matrix, 2) /= size(map_to_cluster)) then
      error stop "cvt_utils: assign_points_to_cluster: size(distance_matrix, 2) /= size(map_to_cluster)."
    end if

    !$OMP parallel do shared(distance_matrix)
    do i_point = 1, size(distance_matrix, 2)
      map_to_cluster(i_point) = minloc(distance_matrix(:, i_point), dim=1)
    end do
    !$OMP end parallel do
  end subroutine assign_points_to_cluster  
end module cvt_utils  