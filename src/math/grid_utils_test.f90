!> Module for unit tests for for the functions in grid_utils.
module grid_utils_test
   use precision, only: dp
   use constants, only: zone
   use modmpi, only: mpiinfo
   use unit_test_framework, only: unit_test_type
   use math_utils, only: all_close
   use grid_utils, only: linspace, grid_3d, phase, fft_frequencies

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

      !> test object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 18

      call test_report%init(n_assertions, mpiglobal)

      ! Run unit tests
      call test_linspace(test_report)

      call test_grid_3d(test_report)

      call test_phase_3d(test_report)

      call test_fft_frequencies(test_report)

      if (present(kill_on_failure)) then
         call test_report%report('grid_utils', kill_on_failure)
      else
         call test_report%report('grid_utils')
      end if

      call test_report%finalise()
   end subroutine grid_utils_test_driver


   !> Test linspace for different inputs.
   subroutine test_linspace(test_report)
      !> Test object
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
      grid = linspace(2.0_dp, 3.0_dp, 5, endpoint=.false.)
      
      call test_report%assert(size(grid) == 5, message='Expect size(grid) == 5')

      call test_report%assert(grid(5) /= 3.0_dp, message='endpoint should not be included')
      
      call test_report%assert(all_close(grid, [2.0_dp, 2.2_dp, 2.4_dp, 2.6_dp, 2.8_dp]), &
                              message='Grid from 2 to 3, containing 5 points but not including the end point')

   end subroutine test_linspace


   !> Test grid_3d for a three dimensional lattice with different numbers of points per dimension.
   !> Expected output: Array of size (2,3,4) with qubicly distributed grid points.
   subroutine test_grid_3d(test_report)
      !> Test object
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
      
                  ref_2d(3, 6) = reshape( [ 1.0_dp, 2.0_dp, 3.0_dp, &   
                                            1.0_dp, 2.0_dp, 3.0_dp, &
                                            1.0_dp, 1.5_dp, 3.0_dp, &
                                            1.0_dp, 1.5_dp, 3.0_dp, &
                                            1.0_dp, 1.0_dp, 3.0_dp, &
                                            1.0_dp, 1.0_dp, 3.0_dp], [3, 6])
      
      grid = grid_3d([ 2, 3, 4 ])
      call test_report%assert(all_close(grid(:, 24), [1.0_dp, 1.0_dp, 1.0_dp]), &
                      'Test grid_3d for [2, 3, 4] with endpoint=.true. does include endpoint.')

      grid = grid_3d([ 2, 3, 4 ], endpoint=.false.)
      call test_report%assert(all_close(grid, ref_3d), &
                      'Test grid_3d for [2, 3, 4] (without end point) as numbers of grid points per dimension. &
                      Expected: 3 x 24 element array that contains the coordinates of the grid points &
                      in a cube with length 1 and [2, 3, 4] number of grid points per dimension.')

      grid = grid_3d([2, 3, 1], [1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 1.0_dp, 1.0_dp])
      call test_report%assert(all_close(grid, ref_2d), &
                      'Test grid_3d for [2, 3, 1] (with end point) as numbers of grid points per dimension. &
                      Expected: 3 x 6 element array that contains the coordinates of the grid points &
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

    phase_grid = phase(grid_3d([2, 2, 1], endpoint = .false.), [0.0_dp, 0.0_dp, 0.0_dp])
    call test_report%assert(all_close(phase_grid, ref_Gamma), &
                    'Test phase_3d at the Gamma point. &
                    Expected: Four element complex vector with complex phase for each grid point at the Gamma point (k=(0, 0, 0)) &
                    for a grid of [2, 2, 1] points per dimension in a cube of lenght 1. &
                    All phases must be (1.0, 0.0).')

    phase_grid = phase(grid_3d([3, 2, 4], endpoint = .false.), [ 0.0_dp, 0.4_dp, 0.6_dp ])
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
&                     Expected output: 3 x 24 element array with the fast fourier freqencies corresponding to the grid in correct order.')
   end subroutine test_fft_frequencies

end module grid_utils_test
