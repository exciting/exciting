module isdf_utils
  use precision, only: dp
  use modmpi, only: mpiinfo
  use xlapack, only: invert_LU, matrix_multiply
  use asserts, only: assert

  implicit none

  private
  public :: isdf

  interface isdf 
    procedure :: isdf_kkp, isdf_samek
  end interface

  contains

  subroutine isdf_kkp(mpi_env, u_reshaped, u_on_interpolation_grid_reshaped, interpolation_point_indices, interpolation_vectors)
    type(mpiinfo), intent(inout) :: mpi_env
    complex(dp), contiguous, intent(in) :: u_reshaped(:, :), u_on_interpolation_grid_reshaped(:, :)
    integer, contiguous, intent(in) :: interpolation_point_indices(:)
    complex(dp), allocatable, intent(out) :: interpolation_vectors(:, :)

    integer :: n_r, n_interpolation_points
    complex(dp), allocatable :: CC_dag_inv(:, :), ZC_dag(:, :), ZC_dag_work(:, :)                      

    call assert(size(u_reshaped, 2) == size(u_on_interpolation_grid_reshaped, 2), &
            'u_reshaped and u_on_interpolation_grid_reshaped have not the same number of combinations: &
                    size(u_reshaped, 2) /= size(u_on_interpolation_grid_reshaped, 2).')

    call assert(size(u_reshaped, 1) >= size(u_on_interpolation_grid_reshaped, 1), &
            'u_reshaped is given on less r-points than u_on_interpolation_grid_reshaped: &
                    size(u_reshaped, 1) < size(u_on_interpolation_grid_reshaped, 1).')
                
    n_r = size(u_reshaped, 1)
    n_interpolation_points = size(u_on_interpolation_grid_reshaped, 1)

    allocate(ZC_dag_work(n_r, n_interpolation_points))

    call matrix_multiply(conjg(u_reshaped), u_on_interpolation_grid_reshaped, ZC_dag_work, trans_B='T')

    ZC_dag = real(conjg(ZC_dag_work) * ZC_dag_work)
    CC_dag_inv = ZC_dag(interpolation_point_indices, :)
    call invert_LU(CC_dag_inv)

    allocate(interpolation_vectors(n_r, n_interpolation_points))
    call matrix_multiply(ZC_dag, CC_dag_inv, interpolation_vectors)

    deallocate(CC_dag_inv, ZC_dag, ZC_dag_work)
  end subroutine isdf_kkp

  subroutine isdf_samek(& 
          mpi_env, &
          u_1_reshaped, &
          u_2_reshaped, &
          u_1_on_interpolation_grid_reshaped, &
          u_2_on_interpolation_grid_reshaped, &
          interpolation_point_indices, &
          interpolation_vectors &
        )
    type(mpiinfo), intent(inout) :: mpi_env
    complex(dp), contiguous :: u_1_reshaped(:, :), u_1_on_interpolation_grid_reshaped(:, :), &
                               u_2_reshaped(:, :), u_2_on_interpolation_grid_reshaped(:, :)
    integer, contiguous, intent(in) :: interpolation_point_indices(:)                           
    complex(dp), allocatable, intent(out) :: interpolation_vectors(:, :)

    integer :: n_r, n_interpolation_points
    complex(dp), allocatable :: ZC_dag(:, :), CC_dag_inv(:, :), ZC_dag_work_1(:, :), ZC_dag_work_2(:, :)

    call assert(size(u_1_reshaped, 1) == size(u_2_reshaped, 1), &
            'u_1_reshaped and u_2_reshaped are not given on the same r-grid: &
                    size(u_1_reshaped, 1) /= size(u_2_reshaped, 1).')

    call assert(size(u_1_on_interpolation_grid_reshaped, 1) == size(u_2_on_interpolation_grid_reshaped, 1), &
            'u_1_on_interpolation_grid_reshaped and u_2_on_interpolation_grid_reshaped are not given on the same r-grid: &
                    size(u_1_on_interpolation_grid_reshaped, 1) /= size(u_2_on_interpolation_grid_reshaped, 1).')

    call assert(size(u_1_reshaped, 1) >= size(u_1_on_interpolation_grid_reshaped, 1), &
            'u_reshaped is given on less r-points than u_on_interpolation_grid_reshaped: &
                    size(u_1_reshaped, 1) < size(u_1_on_interpolation_grid_reshaped, 1).')
                               
    n_r = size(u_1_reshaped, 1)
    n_interpolation_points = size(u_1_on_interpolation_grid_reshaped, 1)

    allocate(ZC_dag_work_1(n_r, n_interpolation_points))
    allocate(ZC_dag_work_2(n_r, n_interpolation_points))

    call matrix_multiply(conjg(u_1_reshaped), u_1_on_interpolation_grid_reshaped, ZC_dag_work_1, trans_B='T')
    call matrix_multiply(conjg(u_2_reshaped), u_2_on_interpolation_grid_reshaped, ZC_dag_work_2, trans_B='T')

    ZC_dag =  conjg(ZC_dag_work_1) * ZC_dag_work_2

    CC_dag_inv = ZC_dag(interpolation_point_indices, :)

    call invert_LU(CC_dag_inv)

    allocate(interpolation_vectors(n_r, n_interpolation_points))
    call matrix_multiply(ZC_dag, CC_dag_inv, interpolation_vectors)

    deallocate(ZC_dag, CC_dag_inv, ZC_dag_work_1, ZC_dag_work_2)
  end subroutine isdf_samek 

end  module isdf_utils

