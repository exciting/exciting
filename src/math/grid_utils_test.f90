!> Module for unit tests for for the functions in grid_utils.
module grid_utils_test
   use precision, only: dp
   use constants, only: zone
   use modmpi, only: mpiinfo
   use unit_test_framework, only: unit_test_type
   use math_utils, only: all_close, all_zero, mod1
   use grid_utils, only: mesh_1d, linspace, concatenate, grid_3d, phase, fft_frequencies, &
                         partial_grid, n_grid_diff, flattened_map, point_in_triangle
    use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices
   implicit none

   private
   public :: grid_utils_test_driver

contains

   !> Run tests for grid utils
   subroutine grid_utils_test_driver(mpiglobal, kill_on_failure)
      !> mpi environment
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program upon failure of an assertion
      logical, intent(in), optional :: kill_on_failure

      !> Test report object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 36

      call test_report%init(n_assertions, mpiglobal)

      ! Run unit tests
      call test_mesh_1d(test_report)

      call test_linspace(test_report)

      call test_concatenate(test_report)

      call test_grid_3d(test_report)

      call test_phase_3d(test_report)

      call test_fft_frequencies(test_report)

      call test_partial_grid(test_report)

      call test_n_grid_diff(test_report)

      call test_flattened_map(test_report)

      if (present(kill_on_failure)) then
         call test_report%report('grid_utils', kill_on_failure)
      else
         call test_report%report('grid_utils')
      end if

      call test_report%finalise()
   end subroutine grid_utils_test_driver


    subroutine test_point_in_triangle(test_report)

      !> Our test object
      type(unit_test_type), intent(inout) :: test_report

      !> First 2d reciprocal lattice vector
      real(dp) :: lat_vec1(2)
      !>  Second 2d reciprocal lattice vector
      real(dp) :: lat_vec2(2)
      !>  High-symmetry point K
      real(dp) :: triangle_center(2)
      !>  High-symmetry point K_bar
      real(dp) :: p_out_1(2)
      !>  Point in the region around the K-valley
      real(dp) :: p_in_1(2)
      !>  Point in the region around the K_bar-valley
      real(dp) :: p_out_2(2)

      lat_vec1 = (/1.0544585551_dp, 0.6087919333_dp/)
      lat_vec2 = (/0._dp, 1.2175838665_dp/)

      triangle_center = 1._dp/3._dp*(lat_vec1 + lat_vec2)
      p_out_1 = 2._dp/3._dp*(lat_vec1 + lat_vec2)

      ! A point that should be in the triangle
      p_in_1 = (/0.2_dp*(lat_vec1(1) +lat_vec2(1)),&
        &  0.16_dp*(lat_vec1(2) + lat_vec2(2))/)
      
      ! A point that should not be in the triangle
      p_out_2 = (/0.75_dp*(lat_vec1(1) + lat_vec2(1)),&
        &  0.5_dp*(lat_vec1(2) + lat_vec2(2))/)

      call test_report%assert(point_in_triangle(lat_vec1, lat_vec2, triangle_center), &
                    & 'Expected:  triangle center is contained in the triangle.')

      call test_report%assert(.not. point_in_triangle(lat_vec1, lat_vec2,&
                          &p_out_1), 'Expected: point not in triangle')

      call test_report%assert(point_in_triangle(lat_vec1, lat_vec2, p_in_1), &
                    & 'Expected: Point in triangle')

      call test_report%assert(.not. point_in_triangle(lat_vec1, lat_vec2,&
                         &p_out_2), 'Expected:  point not in triangle')

    end subroutine

    !> Test mesh_1d
    subroutine test_mesh_1d(test_report)
      !> Unit test object
      type(unit_test_type), intent(inout) :: test_report

      integer, allocatable :: test_range(:), reference_range(:)

      test_range = mesh_1d(4, 4)
      reference_range = [4]
      call test_report%assert(all(test_range == reference_range), 'mesh_1d(4, 4) /= 4.')

      test_range = mesh_1d(1, 5)
      reference_range = [1, 2, 3, 4, 5]
      call test_report%assert(all(test_range == reference_range), &
                            'Test mesh_1d for N_coords as input.')

      test_range = mesh_1d(3, 8)
      reference_range = [3, 4, 5, 6, 7, 8]
      call test_report%assert(all(test_range == reference_range), &
                            'Test mesh_1d for ascending order.')

      test_range = mesh_1d(12, 9)
      reference_range = [12, 11, 10, 9]
      call test_report%assert(all(test_range == reference_range), &
                            'Test mesh_1d for descending order.')

      test_range = mesh_1d(13, 20, 3)
      reference_range = [13, 16, 19]
      call test_report%assert(all(test_range == reference_range), &
                            'Test mesh_1d for ascending order with spacing = 3.')

      test_range = mesh_1d(6, 2, 2)
      reference_range = [6, 4, 2]
      call test_report%assert(all(test_range == reference_range), &
                            'Test mesh_1d for descending order with spacing = 2.')
    end subroutine test_mesh_1d


   !> Test linspace for different inputs.
   subroutine test_linspace(test_report)
      !> Test report object
      type(unit_test_type), intent(inout) :: test_report
      !> Grids
      real(dp), allocatable :: grid(:), descending_grid(:), single_point(:)

      ! Ascending regular grid
      grid = linspace(2.0_dp, 3.0_dp, 0.3_dp)

      call test_report%assert(size(grid) == 4, message='Expect size(grid) == 4')

      call test_report%assert(all_close(grid, [2.0_dp, 2.3_dp, 2.6_dp, 2.9_dp]), &
                              message='Expect grid = [2.0_dp , 2.3_dp, 2.6_dp, 2.9_dp]')

      deallocate (grid)

      ! Descending regular grid
      descending_grid = linspace(3.0_dp, 2.0_dp, 0.3_dp)

      call test_report%assert(size(descending_grid) == 4, &
                              message='Expect size(descending_grid) == 4')

      call test_report%assert(all_close(descending_grid, [3.0_dp, 2.7_dp, 2.4_dp, 2.1_dp]), &
                              message='Expect descending_grid = [3.0_dp , 2.7_dp, 2.4_dp, 2.1_dp]')

      deallocate (descending_grid)

      ! Spacing > Grid limits
      ! This behaviour gets caught in the routine
      ! call test_report%assert(all_close(linsplinspace(2.0_dp, 3.0_dp, 5)ace(2.0_dp, 3.0_dp, 5.3_dp), &
      !                 [2.0_dp]), &
      !                 'Test linspace with spacing larger then the interval. &
      !                 Expected: Array with start as single element ([2.0]).')

      ! Limiting case start = end
      single_point = [2._dp]

      call test_report%assert(all_close(linspace(2.0_dp, 2.0_dp, 5._dp), single_point), &
                              message='If start = end, expect linspace to return start, regardless of spacing')

      deallocate (single_point)

      ! Define number of points in grid, instead of spacing
      grid = linspace(2.0_dp, 3.0_dp, 5)
      
      call test_report%assert(size(grid) == 5, message='Expect size(grid) == 5')

      call test_report%assert(grid(5) == 3._dp, message='endpoint should be included by default')

      call test_report%assert(all_close(grid, [2.0_dp, 2.25_dp, 2.5_dp, 2.75_dp, 3.0_dp]), &
                              message='Grid from 2 to 3, containing 5 points. endpoint included by default')

      deallocate (grid)

      ! Grid defined according to number of points, not including the end point
      grid = linspace(2.0_dp, 3.0_dp, 5, include_last=.false.)
      
      call test_report%assert(size(grid) == 5, message='Expect size(grid) == 5')

      call test_report%assert(grid(5) /= 3.0_dp, message='endpoint should not be included')
      
      call test_report%assert(all_close(grid, [2.0_dp, 2.2_dp, 2.4_dp, 2.6_dp, 2.8_dp]), &
                              message='Grid from 2 to 3, containing 5 points but not including the end point')

   end subroutine test_linspace


   !> Test mesh_1d
   subroutine test_concatenate(test_report)
     !> Unit test object
     type(unit_test_type), intent(inout) :: test_report

     integer, allocatable :: v1_int(:), v2_int(:), v_ref_int(:)
     real(dp), allocatable :: v1_r8(:), v2_r8(:), v_ref_r8(:)
     complex(dp), allocatable :: v1_c8(:), v2_c8(:), v_ref_c8(:)

       v1_int = [1, 2, 4, 3]
       v2_int = [5, 0, 1]
       v_ref_int = [1, 2, 4, 3, 5, 0, 1]
       call test_report%assert(all(concatenate(v1_int, v2_int) ==  v_ref_int), &
               'Test concatonate for integer input.')

       v1_r8 = [1._dp, 2._dp, 0.4_dp, 3._dp, 123._dp]
       v2_r8 = [5._dp]
       v_ref_r8 = [1._dp, 2._dp, 0.4_dp, 3._dp, 123._dp, 5._dp]
       call test_report%assert(all_close(concatenate(v1_r8, v2_r8), v_ref_r8), &
               'Test concatonate for double input.')

       v1_c8 = [cmplx(1., 9., kind=dp), cmplx(-1., 0., kind=dp)]
       v2_c8 = [cmplx(1., 0., kind=dp), cmplx(1., -0.5, kind=dp)]
       v_ref_c8 = [cmplx(1., 9., kind=dp), cmplx(-1., 0., kind=dp), cmplx(1., 0., kind=dp), cmplx(1., -0.5, kind=dp)]
       call test_report%assert(all_close(concatenate(v1_c8, v2_c8), v_ref_c8), &
               'Test concatonate for double complex &input.')
   end subroutine test_concatenate


   !> Test grid_3d for a three dimensional lattice with different numbers of coords per dimension.
   !> Expected output: Array of size (2,3,4) with qubicly distributed grid coords.
   subroutine test_grid_3d(test_report)
      !> Test report object
      type(unit_test_type), intent(inout) :: test_report
      !> test grid
      real(dp), allocatable :: grid(:,:)
      !> reference
      real(dp) :: ref_3d(3, 24) = reshape( [ 0.0_dp,  0.0_dp,              0.0_dp, &
                                             0.5_dp,  0.0_dp,              0.0_dp, &
                                             0.0_dp,  0.3333333333333_dp,  0.0_dp, &
                                             0.5_dp,  0.3333333333333_dp,  0.0_dp, &
                                             0.0_dp,  0.6666666666667_dp,  0.0_dp, &
                                             0.5_dp,  0.6666666666667_dp,  0.0_dp, &
                                             0.0_dp,  0.0_dp,              0.25_dp, &
                                             0.5_dp,  0.0_dp,              0.25_dp, &
                                             0.0_dp,  0.3333333333333_dp,  0.25_dp, &
                                             0.5_dp,  0.3333333333333_dp,  0.25_dp, &
                                             0.0_dp,  0.6666666666667_dp,  0.25_dp, &
                                             0.5_dp,  0.6666666666667_dp,  0.25_dp, &
                                             0.0_dp,  0.0_dp,              0.5_dp, &
                                             0.5_dp,  0.0_dp,              0.5_dp, &
                                             0.0_dp,  0.3333333333333_dp,  0.5_dp, &
                                             0.5_dp,  0.3333333333333_dp,  0.5_dp, &
                                             0.0_dp,  0.6666666666667_dp,  0.5_dp, &
                                             0.5_dp,  0.6666666666667_dp,  0.5_dp, &
                                             0.0_dp,  0.0_dp,              0.75_dp, &
                                             0.5_dp,  0.0_dp,              0.75_dp, &
                                             0.0_dp,  0.3333333333333_dp,  0.75_dp, &
                                             0.5_dp,  0.3333333333333_dp,  0.75_dp, &
                                             0.0_dp,  0.6666666666667_dp,  0.75_dp, &
                                             0.5_dp,  0.6666666666667_dp,  0.75_dp ], [3, 24] ), &
      
                  ref_1d(3, 3) = reshape( [ 1.0_dp, 2.0_dp, 3.0_dp, &
                                            1.0_dp, 1.5_dp, 3.0_dp, &
                                            1.0_dp, 1.0_dp, 3.0_dp ], [3, 3])
      
      grid = grid_3d([ 2, 3, 4 ])
      call test_report%assert(all_close(grid(:, 24), [1.0_dp, 1.0_dp, 1.0_dp]), &
                      'Test grid_3d for [2, 3, 4] with endpoint=.true. does include endpoint.')

      grid = grid_3d([ 2, 3, 4 ], include_last=.false.)
      call test_report%assert(all_close(grid, ref_3d), &
                      'Test grid_3d for [2, 3, 4] (without end point) as numbers of grid points per dimension. &
                      Expected: 3 x 24 element array that contains the coordinates of the grid points &
                      in a cube with length 1 and [2, 3, 4] number of grid points per dimension.')

      grid = grid_3d([1, 3, 1], [1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 1.0_dp, 1.0_dp])
      call test_report%assert(all_close(grid, ref_1d), &
                      'Test grid_3d for [1, 3, 1] (with end point) as numbers of grid points per dimension. &
                      Expected: 3 x 3 element array that contains the coordinates of the grid points &
                      in a cuboid with start = [1.0, 2.0, 3.0] and end = [1.0, 1.0, 1.0].')


    end subroutine test_grid_3d

    !> Test phase_3d at the Gamma point and an arbitrary k-point.
  subroutine test_phase_3d(test_report)
    !> Unit test report
    type(unit_test_type), intent(inout) :: test_report
    !> test phases
    complex(dp) :: phase_single_point
    complex(dp), allocatable :: phase_grid(:)
    !> references
    complex(dp) :: ref_Gamma(4) = [ cmplx( 1.0_dp, 0.0_dp, kind=dp), &
                                    cmplx( 1.0_dp, 0.0_dp, kind=dp), &
                                    cmplx( 1.0_dp, 0.0_dp, kind=dp), &
                                    cmplx( 1.0_dp, 0.0_dp, kind=dp) ] 
      
    complex(dp) :: ref_arb(24) = [ cmplx(  1.0_dp,                   0.0_dp, kind=dp), &
                                   cmplx(  1.0_dp,                   0.0_dp, kind=dp), &
                                   cmplx(  1.0_dp,                   0.0_dp, kind=dp), &
                                   cmplx(  0.30901699437494745_dp, - 0.9510565162951535_dp, kind=dp), &
                                   cmplx(  0.30901699437494745_dp, - 0.9510565162951535_dp, kind=dp), &
                                   cmplx(  0.30901699437494745_dp, - 0.9510565162951535_dp, kind=dp), &
                                   cmplx(  0.5877852522924731_dp,  - 0.8090169943749475_dp, kind=dp), &
                                   cmplx(  0.5877852522924731_dp,  - 0.8090169943749475_dp, kind=dp), &
                                   cmplx(  0.5877852522924731_dp,  - 0.8090169943749475_dp, kind=dp), &
                                   cmplx( -0.5877852522924731_dp,  - 0.8090169943749475_dp, kind=dp), &
                                   cmplx( -0.5877852522924731_dp,  - 0.8090169943749475_dp, kind=dp), &
                                   cmplx( -0.5877852522924731_dp,  - 0.8090169943749475_dp, kind=dp), &
                                   cmplx( -0.30901699437494734_dp, - 0.9510565162951536_dp, kind=dp), &
                                   cmplx( -0.30901699437494734_dp, - 0.9510565162951536_dp, kind=dp), &
                                   cmplx( -0.30901699437494734_dp, - 0.9510565162951536_dp, kind=dp), &
                                   cmplx( -1.0_dp,                 - 1.6653345369377348E-16_dp, kind=dp), &
                                   cmplx( -1.0_dp,                 - 1.6653345369377348E-16_dp, kind=dp), &
                                   cmplx( -1.0_dp,                 - 1.6653345369377348E-16_dp, kind=dp), &
                                   cmplx( -0.9510565162951535_dp,  - 0.3090169943749475_dp, kind=dp), &
                                   cmplx( -0.9510565162951535_dp,  - 0.3090169943749475_dp, kind=dp), &
                                   cmplx( -0.9510565162951535_dp,  - 0.3090169943749475_dp, kind=dp), &
                                   cmplx( -0.5877852522924731_dp,   0.8090169943749473_dp, kind=dp), &
                                   cmplx( -0.5877852522924731_dp,   0.8090169943749473_dp, kind=dp), &
                                   cmplx( -0.5877852522924731_dp,   0.8090169943749473_dp, kind=dp) ]
    
    phase_single_point = phase([0.25_dp, 0.0_dp, 0.0_dp], [1.0_dp, 3.0_dp, 6.0_dp])
    call test_report%assert(all_close(phase_single_point, cmplx(0.0_dp, -1.0_dp, kind=dp)), &
                    'Test phase for single r-point. &
                    Expected: The result is equivalent to exp(- i * pi / 2) = (0.0, -1.0)')

    phase_grid = phase(grid_3d([2, 2, 1], include_last = .false.), [0.0_dp, 0.0_dp, 0.0_dp])
    call test_report%assert(all_close(phase_grid, ref_Gamma), &
                    'Test phase_3d at the Gamma point. &
                    Expected: Four element complex vector with complex phase for each grid point at the Gamma point (k=(0, 0, 0)) &
                    for a grid of [2, 2, 1] points per dimension in a cube of lenght 1. &
                    All phases must be (1.0, 0.0).')

    phase_grid = phase(grid_3d([3, 2, 4], include_last = .false.), [ 0.0_dp, 0.4_dp, 0.6_dp ])
    call test_report%assert(all_close(phase_grid, ref_arb), &
                   'Test phase_3d at an arbitrary k-point. &
                    Expected: 24 element complex vector with complex phase for each grid point at k=(0.0, 0.4, 0.6) &
                    for a grid of [2, 3, 4] points per dimension in a cube of lenght 1.')

  end subroutine test_phase_3d

   !> Test fftfreq for a 3D grid defined by [ 4, 2, 3 ].
   subroutine test_fft_frequencies(test_report)
      !> Unit test report
      type(unit_test_type), intent(inout) :: test_report
      !> reference
      real(dp) :: ref(3, 24) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                        1.0_dp, 0.0_dp, 0.0_dp, &
                                        -2.0_dp, 0.0_dp, 0.0_dp, &
                                        -1.0_dp, 0.0_dp, 0.0_dp, &
                                        0.0_dp, -1.0_dp, 0.0_dp, &
                                        1.0_dp, -1.0_dp, 0.0_dp, &
                                        -2.0_dp, -1.0_dp, 0.0_dp, &
                                        -1.0_dp, -1.0_dp, 0.0_dp, &
                                        0.0_dp, 0.0_dp, 1.0_dp, &
                                        1.0_dp, 0.0_dp, 1.0_dp, &
                                        -2.0_dp, 0.0_dp, 1.0_dp, &
                                        -1.0_dp, 0.0_dp, 1.0_dp, &
                                        0.0_dp, -1.0_dp, 1.0_dp, &
                                        1.0_dp, -1.0_dp, 1.0_dp, &
                                        -2.0_dp, -1.0_dp, 1.0_dp, &
                                        -1.0_dp, -1.0_dp, 1.0_dp, &
                                        0.0_dp, 0.0_dp, -1.0_dp, &
                                        1.0_dp, 0.0_dp, -1.0_dp, &
                                        -2.0_dp, 0.0_dp, -1.0_dp, &
                                        -1.0_dp, 0.0_dp, -1.0_dp, &
                                        0.0_dp, -1.0_dp, -1.0_dp, &
                                        1.0_dp, -1.0_dp, -1.0_dp, &
                                        -2.0_dp, -1.0_dp, -1.0_dp, &
                                        -1.0_dp, -1.0_dp, -1.0_dp], [3, 24])

      call test_report%assert(all_close(fft_frequencies([4, 2, 3]), ref), &
                      'Test fft_frequencies for a 3D grid with [4, 2, 3] numbers of points per dimension. &
                      Expected output: 3 x 24 element array with the fast fourier freqencies corresponding to the grid &
                      in correct order.')
  end subroutine test_fft_frequencies

  !> Test [[n_grid_diff]]
  subroutine test_n_grid_diff(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    integer :: N_ks(3)

    N_ks = [2, 5 ,3]

    call test_report%assert(all(n_grid_diff(N_ks) == [4, 10, 6]), &
                           'Test if n_grid_diff returns the correct n_grid_diff for 3d k-grid. &
                           Expected: [4, 10, 6]')


    N_ks = [4, 1, 23]

    call test_report%assert(all(n_grid_diff(N_ks) == [8, 1, 46]), &
                           'Test if n_grid_diff returns the correct n_grid_diff for 2d k-grid. &
                           Expected: [8, 1, 46]')


    N_ks = [1, 4, 1]

    call test_report%assert(all(n_grid_diff(N_ks) == [1, 8, 1]), &
                           'Test if n_grid_diff returns the correct n_grid_diff for 1d k-grid. &
                           Expected: [1, 8, 1]')


    N_ks = [1, 1, 1]

    call test_report%assert(all(n_grid_diff(N_ks) == [1, 1, 1]), &
                           'Test if n_grid_diff returns the correct n_grid_diff for 10d k-grid. &
                           Expected: [1, 1, 1]')

  end subroutine test_n_grid_diff

  !> Test [[partial_grid]]
  subroutine test_partial_grid(test_report)
    !> Our test object
    type(unit_test_type), intent(inout) :: test_report

    integer, parameter :: N_grid_dense(3) = [4, 4, 1], &
                          N_grid(3) = [2, 2, 1], &
                          int_offset(3) = [0, 0, 0], &
                          indices_3d(3, 4) = reshape([1, 1, 1, &
                                                      2, 1, 1, &
                                                      1, 2, 1, &
                                                      2, 2, 1], [3, 4]), &
                          int_coords_dense(3, 16) = reshape([1, 1, 1, &
                                                             2, 1, 1, &
                                                             3, 1, 1, &
                                                             4, 1, 1, &
                                                             1, 2, 1, &
                                                             2, 2, 1, &
                                                             3, 2, 1, &
                                                             4, 2, 1, &
                                                             1, 3, 1, &
                                                             2, 3, 1, &
                                                             3, 3, 1, &
                                                             4, 3, 1, &
                                                             1, 4, 1, &
                                                             2, 4, 1, &
                                                             3, 4, 1, &
                                                             4, 4, 1], [3, 16])
    integer, allocatable :: map_to_dense_grid(:)
    integer :: divider(3), i, index3d(3)

    divider = N_grid_dense / N_grid
    map_to_dense_grid = partial_grid(N_grid_dense, N_grid)

    call test_report%assert(size(map_to_dense_grid) == product(N_grid), &
            'Test if the map that is given out by partial_grid has the correct size (product(N_grid)).')

    call test_report%assert(all(int_coords_dense(:, map_to_dense_grid) == indices_3d * spread(divider, 2, product(N_grid))), &
            'Test partial grid. Expected: partial_grid gives back a map that maps the points of the big grid to &
            the point of the small grid.')
  end subroutine test_partial_grid

  !> Test flattened_map
  subroutine test_flattened_map(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report

    integer :: ind, i, j, k, N_A_2(2), N_A_3(3), N_B_2(2), N_B_3(3)
    integer, allocatable :: map(:), A_2(:, :), B_2(:, :), A_3(:, :, :), B_3(:, :, :), A_flat(:), B_flat(:)

    N_B_2 = [3, 4]
    allocate(B_2(N_B_2(1), N_B_2(2)))
    do ind=1, product(N_B_2)
      call composite_index_to_indices(ind, N_B_2, i, j)
      B_2(i, j) = ind
    end do

    N_A_2 = [5, 6]
    allocate(A_2(N_A_2(1), N_A_2(2)), source = 0)
    A_2(1 : N_B_2(1), 1 : N_B_2(2)) = B_2


    allocate(A_flat(product(N_A_2)), source = 0)
    B_flat = reshape(B_2, [product(N_B_2)])
    
    map = flattened_map(N_A_2, N_B_2)
    A_flat(map) = B_flat

    call test_report%assert(all(A_2 == reshape(A_flat, N_A_2)), 'Test flattened_map for rank-2 arrays. Expeted: &
            all(A == reshape(A_flat, N_A)).')
    deallocate(A_flat)

    N_B_3 = [13, 7, 153]
    allocate(B_3(N_B_3(1), N_B_3(2), N_B_3(3)))
    do ind=1, product(N_B_3)
      call composite_index_to_indices(ind, N_B_3, i, j, k)
      B_3(i, j, k) = ind
    end do

    N_A_3 = [20, 12, 362]
    allocate(A_3(N_A_3(1), N_A_3(2), N_A_3(3)), source = 0)
    A_3(1 : N_B_3(1), 1 : N_B_3(2), 1 : N_B_3(3)) = B_3

    allocate(A_flat(product(N_A_3)), source = 0)
    B_flat = reshape(B_3, [product(N_B_3)])

    map = flattened_map(N_A_3, N_B_3)
    A_flat(map) = reshape(B_3, [product(N_B_3)])

    call test_report%assert(all(A_3 == reshape(A_flat, N_A_3)), 'Test flattened_map for rank-4 arrays. Expeted: &
            all(A == reshape(A_flat, N_A)).')
  end subroutine test_flattened_map

end module grid_utils_test
