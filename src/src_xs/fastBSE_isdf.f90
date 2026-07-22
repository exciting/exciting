module fastBSE_isdf
  use precision, only: dp, i32
  use math_utils, only: random_order
  use grid_utils, only: mesh_1d, first_element, last_element
  use qrcp_utils, only: qrcp, setup_subsampling_matrix_samek, setup_subsampling_matrix_kkp
  use cvt_utils, only: cvt
  use isdf_utils, only: isdf
  use modmpi, only: mpiinfo, terminate_if_false, terminate_mpi_env
  use modinput, only: input_type

  use xhdf5, only: xhdf5_type, abort_if_not_hdf5
  use xs_hdf5, only: h5ds_wfplot
  use os_utils, only: join_paths
  use xgrid, only: regular_grid_type, setup_unitcell_grid
  use xlapack, only: xgeqp3, qr_column_pivot
  use bse_utils, only: bse_type_to_bool
  use seed_generation, only: set_seed
  use xfftw, only: abort_if_not_fftw3
  use fastBSE_groundstate_properties, only: read_wavefunction_u_hdf5
  use fastBSE_file_strings


  implicit none

  
  private
  public :: fastBSE_isdf_cvt, read_isdf_hdf5 ! , fastBSE_isdf_qrcp


  contains 

  !> Read ISDF coeficients (`[[zeta]]`), the wave functions evaluated on the interpolation points (`[[u_o_isdf]] and/or `[[u_u_isdf]]`).
  !> The coordinates of the real space grid, the wave functions are calculated on (`[[r_vectors]]`) and the indices of the interpolation
  !> points (`[[r_isdf_indices]]`).
  subroutine read_isdf_hdf5(mpi_env, h5file, h5path, zeta, u_o_isdf, u_u_isdf, r_vectors, r_isdf_indices)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Name of the HDF5 file
    type(xhdf5_type), intent(inout) :: h5file
    !> Name of the group in the HDF5 file
    character(*), intent(in) :: h5path
    !> Interpolation coefficients \(\zeta_\mu(\mathbf r)\)
    complex(dp), allocatable, intent(out) :: zeta(:, :)
    !> Occupied wavefunctions, evaluated on the interpolation points.
    complex(dp), allocatable, intent(out), optional :: u_o_isdf(:, :)
    !> Unoccupied wavefunctions, evaluated on the interpolation points.
    complex(dp), allocatable, intent(out), optional :: u_u_isdf(:, :)
    !> Full real space points the wave functions are evaluated on for ISDF
    real(dp), allocatable, intent(out), optional :: r_vectors(:, :)
    !> Indices of the interpolation points \(\mathbf r_\mu\) in the full grid
    integer, allocatable, intent(out), optional :: r_isdf_indices(:)
    
    integer :: n_isdf, n_r, n_ok, n_uk
    integer, allocatable :: shape_r(:), shape_indices(:), shape_zeta(:), shape_u_o_isdf(:), shape_u_u_isdf(:)
    character(:), allocatable :: group 

    group = join_paths(h5path, isdf_group)
    
    ! Read data for the exchange kernel
    if (present(u_o_isdf) .and. present(u_u_isdf)) then
      group = join_paths(group, vexc_ou_group)

    ! Read data for the occupied part of the screened kernel
    else if (present(u_o_isdf) .and. (.not. present(u_u_isdf))) then
      group = join_paths(group, wscr_oo_group)

    ! Read data for the unoccupied part of the screened kernel
    else if ((.not. present(u_o_isdf)) .and. present(u_u_isdf)) then
      group = join_paths(group, wscr_uu_group)

    ! You did something wrong
    else 
      call terminate_mpi_env(mpi_env, 'Error: read_fastBSE_isdf_vexc_hdf5: At least one of u_o_isdf or u_u_isdf must be present.')
    end if 

    call h5file%dataset_shape(group, isdf_indices_dataset, shape_indices)
    call h5file%dataset_shape(group, zeta_dataset, shape_zeta, complex_dataset = .true.)
    call h5file%dataset_shape(group, rspace_coordinates_dataset, shape_r)

    n_isdf = shape_indices(1)
    n_r = shape_r(2)

    call terminate_if_false(mpi_env, shape_zeta(1) == n_r, &
            'read_fastBSE_isdf_vexc_hdf5: For dataset ' // zeta_dataset // ' dimension 1 is not &
            the same as the size of dataset ' // rspace_coordinates_dataset //' (n_r), as expected.')

    call terminate_if_false(mpi_env, shape_zeta(2) == n_isdf, &
            'read_fastBSE_isdf_vexc_hdf5: For dataset ' // zeta_dataset // ' dimension 2 is not &
            the same as the size of dataset ' // isdf_indices_dataset //' (n_isdf), as expected.')

    
    allocate(zeta(n_r, n_isdf))

    call h5file%read(group, zeta_dataset, zeta, [1, 1])

    if (present(u_o_isdf)) then
      call h5file%dataset_shape(group, u_o_isdf_dataset, shape_u_o_isdf, complex_dataset = .true.)
      call terminate_if_false(mpi_env, shape_u_o_isdf(1) == n_isdf, &
              'read_fastBSE_isdf_vexc_hdf5: For dataset ' // u_o_isdf_dataset // ' dimension 1 is not &
              the same as the size of dataset ' // isdf_indices_dataset //' (n_isdf), as expected.')

      n_ok = shape_u_o_isdf(2)
      allocate(u_o_isdf(n_isdf, n_ok))
      call h5file%read(group, u_o_isdf_dataset, u_o_isdf)
    end if

    if (present(u_u_isdf)) then
      call h5file%dataset_shape(group, u_u_isdf_dataset, shape_u_u_isdf, complex_dataset = .true.)
      call terminate_if_false(mpi_env, shape_u_u_isdf(1) == n_isdf, &
              'read_fastBSE_isdf_vexc_hdf5: For dataset ' // u_u_isdf_dataset // ' dimension 1 is not &
              the same as the size of dataset ' // isdf_indices_dataset //' (n_isdf), as expected.')

      n_uk = shape_u_u_isdf(2)
      allocate(u_u_isdf(n_isdf, n_uk))
      call h5file%read(group, u_u_isdf_dataset, u_u_isdf)
    end if

    if (present(r_isdf_indices)) then
      allocate(r_isdf_indices(n_isdf))
      call h5file%read(group, isdf_indices_dataset, r_isdf_indices)
    end if

    if (present(r_vectors)) then
      allocate(r_vectors(3, n_r))
      call h5file%read(group, rspace_coordinates_dataset, r_vectors)
    end if 
  end subroutine 

  !> Calculate ISDF with interpolation points obtained by [[cvt]].
  subroutine fastBSE_isdf_cvt(mpi_env, input, h5file, h5group, info_unit)
    type(mpiinfo), intent(in) :: mpi_env
    type(input_type), intent(in) :: input
    character(*), intent(in) :: h5file, h5group
    integer, intent(in) :: info_unit

    ! HDF5
    type(xhdf5_type) :: h5
    character(:), allocatable :: group_isdf, group

    character(:), allocatable :: seed, bse_type
    logical :: calculate_vexc, calculate_wscr, has_converged

    integer :: n_r, n_ok, n_uk, n_k, n_isdf, r_sampling(3), n_combinations, nisdf(3), i_function
    integer :: cvtsteplim, steps
    real(dp) :: lattice(3, 3), r_offset(3), epslat, time_start, time_end, time_cvt, time_isdf, tolerance, deviation

    
    integer, allocatable :: index_map_u(:,:), index_map_o(:, :)
    complex(dp), allocatable :: u_o(:, :), u_u(:, :), u_o_isdf(:, :), u_u_isdf(:, :)
    
    ! CVT + ISDF
    integer, allocatable :: r_isdf_indices(:)
    real(dp), allocatable :: rho(:), real_space_coordinates(:, :)
    complex(dp), allocatable :: M(:, :), tau(:), zeta(:, :)

    type(regular_grid_type) :: r_grid 

    ! Intialize input data
    seed          = input%xs%fastBSE%seed
    lattice       = input%structure%crystal%basevect
    epslat        = input%structure%epslat
    r_offset      = spread(0._dp, 1, 3)
    r_sampling    = input%xs%fastBSE%ngridr
    r_grid        = setup_unitcell_grid(r_sampling, r_offset, lattice, epslat)
    n_r           = r_grid%number_of_points()
    n_k           = product(input%xs%ngridk)
    cvtsteplim = input%xs%fastBSE%cvtsteplim
    tolerance     = input%xs%fastBSE%cvttol
    nisdf = input%xs%fastBSE%nisdf
    bse_type      = input%xs%BSE%bsetype

    call abort_if_not_fftw3(mpi_env, "Error(fastBSE_isdf_cvt): exciting needs to be linked to FFTW3 for running fastBSE.")
    call abort_if_not_hdf5(mpi_env, "Error(fastBSE_isdf_cvt): exciting needs to be compiled with HDF5 to run fastBSE module.")
    
    call bse_type_to_bool(bse_type, calculate_vexc, calculate_wscr)

    ! Only calculate ISDF if the interaction kernels are needed.
    if (.not. (calculate_vexc .or. calculate_wscr)) return
    
    call set_seed(seed)

    call h5%initialize(h5file, mpi_env)
    call h5%initialize_group(h5group, isdf_group)
    group_isdf = join_paths(h5group, isdf_group)
    call read_wavefunction_u_hdf5(mpi_env, input, h5, h5group, u_u, u_o)

    n_uk = size(u_u, 2)
    n_ok = size(u_o, 2)
    rho = sqrt(sum(abs(u_o) ** 2, dim=2) + sum(abs(u_u) ** 2, dim=2))
    real_space_coordinates = r_grid%coordinate_array()

    if(calculate_vexc) then
      
      ! Calculate interpolation points via CVT
      call timesec(time_start)
      n_isdf = min(nisdf(1), n_ok * n_uk / n_k)
      steps = cvtsteplim
      r_isdf_indices = random_order(n_r, n_out=n_isdf)
      call cvt(real_space_coordinates, rho, tolerance, steps, r_isdf_indices, has_converged, deviation)
      call timesec(time_end)
      time_cvt = time_end - time_start

      ! Calculate ISDF zeta zeta
      call timesec(time_start)
      u_o_isdf = u_o(r_isdf_indices, :)
      u_u_isdf = u_u(r_isdf_indices, :)
      call isdf(mpi_env, u_o, u_u, u_o_isdf, u_u_isdf, r_isdf_indices, zeta)
      call timesec(time_end)
      time_isdf = time_end - time_start

      ! Write results to file
      call h5%initialize_group(group_isdf, vexc_ou_group)
      group = join_paths(group_isdf, vexc_ou_group)
      call h5%write(group, isdf_indices_dataset, r_isdf_indices)
      call h5%write(group, zeta_dataset, zeta)
      call h5%write(group, rspace_coordinates_dataset, r_grid%coordinate_array())
      call h5%write(group, u_o_isdf_dataset, u_o_isdf)
      call h5%write(group, u_u_isdf_dataset, u_u_isdf)

      call write_info(info_unit, "V_exc", time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)

      ! Clean up
      deallocate(r_isdf_indices, zeta, u_o_isdf, u_u_isdf)
    
    end if

    if (calculate_wscr) then 

      ! occupied states

      ! Calculate interpolation points via CVT
      call timesec(time_start)
      n_isdf = min(nisdf(2), n_ok**2)
      steps = cvtsteplim
      r_isdf_indices = random_order(n_r, n_isdf)
      call cvt(real_space_coordinates, rho, tolerance, steps, r_isdf_indices, has_converged, deviation)
      call timesec(time_end)
      time_cvt = time_end - time_start

      ! Calculate ISDF zeta zeta
      call timesec(time_start)
      u_o_isdf = u_o(r_isdf_indices, :)
      call isdf(mpi_env, u_o, u_o_isdf, r_isdf_indices, zeta)
      call timesec(time_end)
      time_isdf = time_end - time_start
      
      ! Write results to file
      call h5%initialize_group(group_isdf, wscr_oo_group)
      group = join_paths(group_isdf, wscr_oo_group)
      call h5%write(group, isdf_indices_dataset, r_isdf_indices)
      call h5%write(group, zeta_dataset, zeta)
      call h5%write(group, rspace_coordinates_dataset, r_grid%coordinate_array())
      call h5%write(group, u_o_isdf_dataset, u_o_isdf)
      
      ! Write Info to file
      call write_info(info_unit, "W_scr occupied", time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)

      ! Clean up
      deallocate(r_isdf_indices, zeta, u_o_isdf)
      
      ! unoccupied states

      ! Calculate interpolation points via CVT
      call timesec(time_start)
      n_isdf = min(nisdf(3), n_uk**2)
      steps = cvtsteplim
      r_isdf_indices = random_order(n_r, n_isdf)
      call cvt(real_space_coordinates, rho, tolerance, steps, r_isdf_indices, has_converged, deviation)
      call timesec(time_end)
      time_cvt = time_end - time_start

      ! Calculate ISDF zeta zeta
      call timesec(time_start)
      u_u_isdf = u_u(r_isdf_indices, :)
      call isdf(mpi_env, u_u, u_u_isdf, r_isdf_indices, zeta)
      call timesec(time_end)
      time_isdf = time_end - time_start

      ! Write results to file
      call h5%initialize_group(group_isdf, wscr_uu_group)
      group = join_paths(group_isdf, wscr_uu_group)
      call h5%write(group, isdf_indices_dataset, r_isdf_indices)
      call h5%write(group, zeta_dataset, zeta)
      call h5%write(group, rspace_coordinates_dataset, r_grid%coordinate_array())
      call h5%write(group, u_u_isdf_dataset, u_u_isdf)

      ! Write Info to file
      call write_info(info_unit, "W_scr unoccupied", time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)

      ! Clean up
      deallocate(r_isdf_indices, zeta, u_u_isdf)
    end if

    call h5%finalize()

    deallocate(u_o, u_u)   
    
    contains

    subroutine write_info(info_unit, taskname, time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)
      integer, intent(in) :: info_unit
      character(*), intent(in) :: taskname
      real(dp) :: time_cvt, time_isdf
      logical :: has_converged 
      integer :: steps
      real(dp) :: deviation
      integer :: n_isdf
      
      write(info_unit, *)
      write(info_unit, "(A, A)")     "ISDF + CVT done for ", taskname  
      
      if(has_converged) Then
        write(info_unit, "(A, I4, A)") "CVT procedure has converged in ", steps, " steps."
      else 
        write(info_unit, "(A, I4, A)") "CVT procedure has not converged in ", steps, " steps."
      end if 
      write(info_unit, "(A, ES21.14)") "Deviation :", deviation
      write(info_unit, *)
      
      write(info_unit, "(A, I8)")    "number of interpolation points: ",  n_isdf
      write(info_unit, *)
  
      write(info_unit, "(A, F15.6)") "Time(CVT)  (s): ", time_cvt
      write(info_unit, "(A, F15.6)") "Time(ISDF) (s): ", time_isdf
      write(info_unit, *)
      
    end subroutine write_info

  end subroutine fastBSE_isdf_cvt

end module fastBSE_isdf