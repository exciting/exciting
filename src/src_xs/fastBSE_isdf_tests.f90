module fastBSE_isdf_tests
  use precision, only: dp 
  use math_utils, only: random_order
  use qrcp_utils, only: qrcp, setup_subsampling_matrix_samek, setup_subsampling_matrix_kkp
  use cvt_utils, only: cvt
  use isdf_utils, only: isdf
  use modmpi, only: mpiinfo
  use modinput, only: input_type
  use asserts, only: assert
  use xhdf5, only: xhdf5_type
  use xs_hdf5, only: h5ds_wfplot
  use os_utils, only: join_paths
  use xgrid, only: regular_grid_type, setup_unitcell_grid
  use xlapack, only: xgeqp3, qr_column_pivot
  use bse_utils, only: bse_type_to_bool
  use seed_generation, only: set_seed
  use isdf_test_utils, only: calculate_wavefunction_product_samek, calculate_isdf_wavefunction_product_samek
  use fastBSE_isdf, only: h5group_isdf_vexc, &
                                    h5group_isdf_wscr_o, &
                                    h5group_isdf_wscr_u

  contains 


  subroutine fastBSE_isdf_vexc_test(mpi_env, input, h5file, h5group, info_unit)
    use fastBSE_isdf, only: h5group_isdf_vexc, h5ds_zeta, &
                                      h5ds_wfplot_isdf_o, &
                                      h5ds_wfplot_isdf_u
    type(mpiinfo), intent(inout) :: mpi_env
    type(input_type), intent(in) :: input
    character(*), intent(in) :: h5file, h5group
    integer, intent(in) :: info_unit
    

    ! HDF5
    type(xhdf5_type) :: h5
    character(:), allocatable :: path

    character(:), allocatable :: seed, bse_type
    logical :: calculate_vexc, calculate_wscr

    integer, allocatable :: dims(:)
    integer :: n_r, n_o, n_u, n_k, n_isdf
    integer :: i 
    real(dp) :: time_start, time_end, error, max_value
    real(dp), allocatable :: error_array(:, :, :, :)
    complex(dp), allocatable :: u_o(:, :, :), u_u(:, :, :), u_o_isdf(:, :, :), u_u_isdf(:, :, :), zeta(:, :), &
                                wf_product(:, :, :, :), wf_product_isdf(:, :, :, :)
    


    call h5%initialize(h5file, mpi_env%comm)

    ! Intialize input data
    n_r           = product(input%xs%fastBSE%rsampling)
    n_o           = input%xs%BSE%nstlbse(2) - input%xs%BSE%nstlbse(1) + 1
    n_u           = input%xs%BSE%nstlbse(4) - input%xs%BSE%nstlbse(3) + 1
    n_k           = product(input%xs%ngridk)
    call h5%dataset_shape(h5group_isdf_vexc, h5ds_zeta, dims, complex_dataset=.true.)
    n_r = dims(1)
    n_isdf = dims(2)

    ! Read wavefunctions
    allocate(u_o(n_r, n_o, n_k))
    allocate(u_u(n_r, n_u, n_k))
    allocate(u_o_isdf(n_isdf, n_o, n_k))
    allocate(u_u_isdf(n_isdf, n_u, n_k))
    allocate(zeta(n_r, n_isdf))
    
    path = h5group
    call h5%read(h5group, h5ds_wfplot, u_o, [1, 1, 1])
    call h5%read(h5group, h5ds_wfplot, u_u, [1, n_o + 1, 1])

    path = join_paths(h5group, h5group_isdf_vexc)
    call h5%read(path, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1, 1])
    call h5%read(path, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1, 1])
    call h5%read(path, h5ds_zeta, zeta, [1, 1, 1])

    ! Setup exact wavefunction products as reference
    call timesec(time_start)
    allocate(wf_product(n_r, n_o, n_u, n_k))
    call calculate_wavefunction_product_samek(u_o, u_u, wf_product)
    call timesec(time_end)
    write(info_unit, *)
    write(info_unit, '(A, F15.6)') 'Setup wavefunction product. Time spent (s) ', time_end - time_start

    ! Setup wavefunction products with ISDF
    call timesec(time_start)
    allocate(wf_product_isdf(n_r, n_o, n_u, n_k))
    call calculate_isdf_wavefunction_product_samek(zeta, u_o_isdf, u_u_isdf, wf_product_isdf)
    call timesec(time_end)
    write(info_unit, *)
    write(info_unit, '(A, F15.6)') 'Setup wavefunction product from ISDF. Time spent (s) ', time_end - time_start

    ! Calculate deveation array
    call timesec(time_start)
    max_value = max(maxval(abs(wf_product)), maxval(abs(wf_product_isdf)))
    error_array = abs(wf_product - wf_product_isdf)
    call timesec(time_end)
    write(info_unit, *)
    write(info_unit, '(A, F15.10)') 'Calculated error array. Time spent (s) ', time_end - time_start
    write(info_unit, '(A, F15.10)') 'Max. error          : ', maxval(error_array)
    write(info_unit, '(A, F15.10)') 'Max. relative error : ', maxval(error_array / max_value)
    write(info_unit, '(A, F15.10)') 'MSE                 : ', sum(error_array**2) / (n_r * n_o * n_u * n_k)
    write(info_unit, '(A, F15.10)') 'Relative MSE        : ', sum((error_array / max_value)**2) / (n_r * n_o * n_u * n_k)
    write(info_unit, *)
    
  end subroutine

end module fastBSE_isdf_tests 
