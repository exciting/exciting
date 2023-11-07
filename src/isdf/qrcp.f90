module qrcp_utils
  use precision, only: dp
  use constants, only: zzero 
  use xfftw, only: fft_type, FFTW_FORWARD
  use asserts, only: assert
  use multi_index_conversion, only: composite_index_to_indices, indices_to_composite_index
  use xlapack, only: outer_product, qr_column_pivot
  use modmpi, only: mpiinfo
  use math_utils, only: random_order, random_units
  use exciting_mpi, only: xmpi_bcast

  implicit none

  private
  public :: qrcp, setup_subsampling_matrix_samek, setup_subsampling_matrix_kkp

  interface qrcp
    procedure :: qrcp_samek, qrcp_kkp
  end interface

  contains

  subroutine qrcp_samek(mpi_env, n_interpolation_points, n_sub, u_1, u_2, interpolation_point_indices)
    type(mpiinfo), intent(inout) :: mpi_env
    integer, intent(in) :: n_sub, n_interpolation_points
    complex(dp), contiguous, intent(in) :: u_1(:, :, :), u_2(:, :, :)
    integer, allocatable, intent(out) :: interpolation_point_indices(:)

    integer :: n_r, n_combinations
    complex(dp), allocatable :: M(:, :)
    integer, allocatable :: r_work(:)

    call assert(size(u_1, 1) == size(u_2, 1), &
            'u_1 and u_2 are not given on the same r-grid: size(u_1, 1) /= size(u_2, 1).')

    call assert(size(u_1, 3) == size(u_2, 3), &
            'u_1 and u_2 are not given on the same k-grid: size(u_1, 3) /= size(u_2, 3).')

    n_r = size(u_1, 1)
    n_combinations = size(u_1, 2) * size(u_2, 2) * size(u_1, 3)

    call assert(n_r >= n_interpolation_points, &
            'The given number of interpolation points is larger than the number of r-points: &
                    n_r < n_interpolation_points.')

    call assert(n_combinations >= n_sub, &
            'The number of subsampling rows is larger than the number of possible combination between u_1 and u_2 with the same k-point: &
                    n_combinations < n_sub.')

    call setup_subsampling_matrix_samek(mpi_env, n_sub, u_1, u_2, M)

    allocate(r_work(n_r))
    call qr_column_pivot(M, r_work)
    interpolation_point_indices = r_work(:n_interpolation_points)

    deallocate(M, r_work)
  end subroutine qrcp_samek


  subroutine qrcp_kkp(mpi_env, n_interpolation_points, n_sub, u, interpolation_point_indices)
    type(mpiinfo), intent(inout) :: mpi_env
    integer, intent(in) :: n_sub, n_interpolation_points
    complex(dp), contiguous, intent(in) :: u(:, :, :)
    integer, allocatable, intent(out) :: interpolation_point_indices(:)

    integer :: n_r, n_combinations
    complex(dp), allocatable :: M(:, :)
    integer, allocatable :: r_work(:)

    n_r = size(u, 1)
    n_combinations = size(u, 2) * size(u, 3) ! n_bands * n_k

    call assert(n_r >= n_interpolation_points, &
            'The given number of interpolation points is larger than the number of r-points: &
                    n_r < n_interpolation_points.')

    call assert(n_combinations >= n_sub, &
            'The number of subsampling rows is larger than the number of possible combination between u and itself: &
                    n_combinations < n_sub.')

    call setup_subsampling_matrix_kkp(mpi_env, n_sub, u, M)

    allocate(r_work(n_r))
    call qr_column_pivot(M, r_work)
    interpolation_point_indices = r_work(:n_interpolation_points)
    
    deallocate(M, r_work)
  end subroutine qrcp_kkp


  subroutine setup_subsampling_matrix_samek(mpi_env, n_sub, u_1, u_2, M)
    type(mpiinfo), intent(inout) :: mpi_env
    integer, intent(in) :: n_sub
    complex(dp), intent(in) :: u_1(:, :, :), u_2(:, :, :)
    complex(dp), allocatable, intent(out) :: M(:, :)

    integer :: n_r, n_bands_1, n_bands_2, n_k, n_combinations, r_index, band_index_1, band_index_2, composite_band_index
    complex(dp), allocatable :: uu_at_r(:, :), uu_at_r_reshaped(:), random_unit(:)
    integer, allocatable :: random_subsampling_order(:)
    type(fft_type):: fft

    n_r = size(u_1, 1)
    n_bands_1 = size(u_1, 2)
    n_bands_2 = size(u_2, 2)
    n_k = size(u_1, 3)
    n_combinations = n_bands_1 * n_bands_2 * n_k

    if(allocated(M)) deallocate(M)
    allocate(M(n_sub, n_r))

    
    if(mpi_env%is_root) then
      random_unit = random_units(n_combinations)
      random_subsampling_order = random_order(n_combinations)
      random_subsampling_order = random_subsampling_order(: n_sub)
    end if

    call xmpi_bcast(mpi_env, random_unit)
    call xmpi_bcast(mpi_env, random_subsampling_order)

    allocate(uu_at_r(n_bands_1 * n_bands_2, n_k))
    allocate(uu_at_r_reshaped(n_combinations))
    call fft%initialize([n_combinations], FFTW_FORWARD, uu_at_r_reshaped)

    do r_index = 1, n_r

      do composite_band_index = 1, n_bands_1 * n_bands_2
        call composite_index_to_indices(composite_band_index, [n_bands_1, n_bands_2], band_index_1, band_index_2)
        uu_at_r(composite_band_index, :) = conjg(u_1(r_index, band_index_1, :)) * u_2(r_index, band_index_2, :)
      end do

      uu_at_r_reshaped = random_unit * reshape(uu_at_r, [n_combinations])

      call fft%execute(uu_at_r_reshaped)

      M(:, r_index) = uu_at_r_reshaped(random_subsampling_order)
    end do

    deallocate(uu_at_r, uu_at_r_reshaped, random_unit, random_subsampling_order)
  end subroutine setup_subsampling_matrix_samek

  subroutine setup_subsampling_matrix_kkp(mpi_env, n_sub, u, M)
    type(mpiinfo), intent(inout) :: mpi_env
    integer, intent(in) :: n_sub
    complex(dp), intent(in) :: u(:, :, :)
    complex(dp), intent(out), allocatable :: M(:, :)

    integer :: n_r, n_bands, n_k, n_combinations, r_index
    integer, allocatable :: random_subsampling_order(:)
    complex(dp), allocatable :: u_at_r(:), uu_at_r(:, :), random_unit(:)
    type(fft_type):: fft

    n_r = size(u, 1)
    n_bands = size(u, 2)
    n_k = size(u, 3)
    n_combinations = n_k * n_bands

    if(allocated(M)) deallocate(M)
    allocate(M(n_sub**2, n_r))

    

    if(mpi_env%is_root) then
      random_unit = random_units(n_combinations)
      random_subsampling_order = random_order(n_combinations)
      random_subsampling_order = random_subsampling_order(: n_sub)
    end if

    call xmpi_bcast(mpi_env, random_unit)
    call xmpi_bcast(mpi_env, random_subsampling_order) 

    allocate(uu_at_r(n_sub, n_sub))
    allocate(u_at_r(n_combinations))
    call fft%initialize([n_combinations], FFTW_FORWARD, u_at_r)

    do r_index = 1, n_r
      u_at_r = random_unit * reshape(u(r_index, :, :), [n_combinations])

      call fft%execute(u_at_r)

      u_at_r = u_at_r(random_subsampling_order)

      uu_at_r = zzero 
      call outer_product(conjg(u_at_r), u_at_r, uu_at_r)

      M(:, r_index) = reshape(uu_at_r, [n_sub**2])
    end do

    deallocate(u_at_r, uu_at_r, random_unit, random_subsampling_order)
  end subroutine setup_subsampling_matrix_kkp

end module qrcp_utils