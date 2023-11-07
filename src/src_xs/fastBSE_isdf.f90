module fastBSE_isdf
  use precision, only: dp, i32
  use math_utils, only: random_order
  use grid_utils, only: mesh_1d, first_element, last_element
  use qrcp_utils, only: qrcp, setup_subsampling_matrix_samek, setup_subsampling_matrix_kkp
  use cvt_utils, only: cvt
  use isdf_utils, only: isdf
  use modmpi, only: mpiinfo
  use modinput, only: input_type
  use asserts, only: assert
  use xhdf5, only: xhdf5_type, abort_if_not_hdf5
  use xs_hdf5, only: h5ds_wfplot
  use os_utils, only: join_paths
  use xgrid, only: regular_grid_type, setup_unitcell_grid
  use xlapack, only: xgeqp3, qr_column_pivot
  use bse_utils, only: bse_type_to_bool
  use seed_generation, only: set_seed
  use xfftw, only: abort_if_not_fftw3
  use fastBSE_write_wfplot, only: h5group_wfplot, h5ds_u, h5ds_r_sampling, h5ds_k_list, h5ds_band_list
  use fastBSE_transitions, only: h5group_transitions, h5ds_uo_limits, h5ds_band_idx


  implicit none

  private
  public :: fastBSE_isdf_cvt! , fastBSE_isdf_qrcp

  character(*), parameter, public :: h5group_isdf_vexc = "isdf_vexc"
  character(*), parameter, public :: h5group_isdf_wscr_o = "isdf_wscr_occupied"
  character(*), parameter, public :: h5group_isdf_wscr_u = "isdf_wscr_unoccupied"

  character(*), parameter, public :: h5ds_points = "grid_points"
  character(*), parameter, public :: h5ds_indices = "indices"
  character(*), parameter, public :: h5ds_zeta = "zeta"
  character(*), parameter, public :: h5ds_wfplot_isdf_o = "wfplot_isdf_occ"
  character(*), parameter, public :: h5ds_wfplot_isdf_u = "wfplot_isdf_unocc"

  contains 






  !> Calculate ISDF with interpolation points obtained by [[]]
  subroutine fastBSE_isdf_cvt(mpi_env, input, h5file, h5group, info_unit)
    type(mpiinfo), intent(inout) :: mpi_env
    type(input_type), intent(in) :: input
    character(*), intent(in) :: h5file, h5group
    integer, intent(in) :: info_unit

    ! HDF5
    type(xhdf5_type) :: h5
    character(:), allocatable :: group

    character(:), allocatable :: seed, bse_type
    logical :: calculate_vexc, calculate_wscr, has_converged

    integer :: n_r, n_ok, n_uk, n_k, n_isdf, r_sampling(3), n_combinations, nisdf(3), i_function
    integer :: cvtsteplim, steps
    real(dp) :: lattice(3, 3), r_offset(3), epslat, time_start, time_end, time_cvt, time_isdf, tolerance, deviation

    
    integer, allocatable :: index_map_u(:,:), index_map_o(:, :)
    complex(dp), allocatable :: u_o(:, :), u_u(:, :), u_o_isdf(:, :), u_u_isdf(:, :)
    
    ! CVT + ISDF
    integer, allocatable :: indices(:)
    real(dp), allocatable :: rho(:), points(:, :)
    complex(dp), allocatable :: M(:, :), tau(:), zeta(:, :)

    type(regular_grid_type) :: r_grid 

    ! Intialize input data
    seed          = input%xs%fastBSE%seed
    lattice       = input%structure%crystal%basevect
    epslat        = input%structure%epslat
    r_offset      = spread(0._dp, 1, 3)
    r_sampling    = input%xs%fastBSE%rsampling
    r_grid        = setup_unitcell_grid(r_sampling, r_offset, lattice, epslat)
    n_r           = r_grid%number_of_points()
    n_k           = product(input%xs%ngridk)
    cvtsteplim = input%xs%fastBSE%cvtsteplim
    tolerance     = input%xs%fastBSE%cvttol
    nisdf = input%xs%fastBSE%nisdf
    bse_type      = input%xs%BSE%bsetype

    call abort_if_not_fftw3(mpi_env, "Error(fastBSE_write_u): exciting needs to be linked to FFTW3 for running fastBSE.")
    call abort_if_not_hdf5(mpi_env, "Error(fastBSE_write_u): exciting needs to be compiled with HDF5 to run fastBSE module.")
    
    call bse_type_to_bool(bse_type, calculate_vexc, calculate_wscr)

    if (.not. (calculate_vexc .or. calculate_wscr)) return
    
    call set_seed(seed)

    call load_u(mpi_env, input, h5file, h5group, u_u, u_o)

    n_uk = size(u_u, 2)
    n_ok = size(u_o, 2)
    
    rho = sqrt(sum(abs(u_o) ** 2, dim=2) + sum(abs(u_u) ** 2, dim=2))
    points = r_grid%coordinate_array()

    call h5%initialize(h5file, mpi_env%comm)

    if(calculate_vexc) then
      
      ! Calculate interpolation points via CVT
      call timesec(time_start)
      n_isdf = min(nisdf(1), n_ok * n_uk / n_k)
      steps = cvtsteplim
      indices = random_order(n_r, n_out=n_isdf)
      call cvt(points, rho, tolerance, steps, indices, has_converged, deviation)
      call timesec(time_end)
      time_cvt = time_end - time_start

      ! Calculate ISDF coefficients zeta
      call timesec(time_start)
      u_o_isdf = u_o(indices, :)
      u_u_isdf = u_u(indices, :)
      call isdf(mpi_env, u_o, u_u, u_o_isdf, u_u_isdf, indices, zeta)
      call timesec(time_end)
      time_isdf = time_end - time_start

      ! Write results to file
      call h5%initialize_group(h5group, h5group_isdf_vexc)
      group = join_paths(h5group, h5group_isdf_vexc)
      
      call h5%write(group, h5ds_indices, indices, [1], [n_isdf])
      call h5%write(group, h5ds_zeta, zeta, [1, 1], [r_grid%number_of_points(), n_isdf])
      call h5%write(group, h5ds_points, r_grid%coordinate_array(), [1, 1], [3, r_grid%number_of_points()])
      call h5%write(group, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1], shape(u_o_isdf))
      call h5%write(group, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1], shape(u_u_isdf))

      call write_info(info_unit, "V_exc", time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)

      ! Clean up
      deallocate(indices, zeta, u_o_isdf, u_u_isdf)
    
    end if

    if (calculate_wscr) then 

      ! occupied states

      ! Calculate interpolation points via CVT
      call timesec(time_start)
      n_isdf = min(nisdf(2), n_ok**2)
      steps = cvtsteplim
      indices = random_order(n_r, n_isdf)
      call cvt(points, rho, tolerance, steps, indices, has_converged, deviation)
      call timesec(time_end)
      time_cvt = time_end - time_start

      ! Calculate ISDF coefficients zeta
      call timesec(time_start)
      u_o_isdf = u_o(indices, :)
      call isdf(mpi_env, u_o, u_o_isdf, indices, zeta)
      call timesec(time_end)
      time_isdf = time_end - time_start
      
      ! Write results to file
      call h5%initialize_group(h5group, h5group_isdf_wscr_o)
      group = join_paths(h5group, h5group_isdf_wscr_o)
      call h5%write(group, h5ds_indices, indices, [1], [n_isdf])
      call h5%write(group, h5ds_zeta, zeta, [1, 1], [r_grid%number_of_points(), n_isdf])
      call h5%write(group, h5ds_points, r_grid%coordinate_array(), [1, 1], [3, r_grid%number_of_points()])
      call h5%write(group, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1], shape(u_o_isdf))
      
      ! Write Info to file
      call write_info(info_unit, "W_scr occupied", time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)

      ! Clean up
      deallocate(indices, zeta, u_o_isdf)
      
      ! unoccupied states

      ! Calculate interpolation points via CVT
      call timesec(time_start)
      n_isdf = min(nisdf(3), n_uk**2)
      steps = cvtsteplim
      indices = random_order(n_r, n_isdf)
      call cvt(points, rho, tolerance, steps, indices, has_converged, deviation)
      call timesec(time_end)
      time_cvt = time_end - time_start

      ! Calculate ISDF coefficients zeta
      call timesec(time_start)
      u_u_isdf = u_u(indices, :)
      call isdf(mpi_env, u_u, u_u_isdf, indices, zeta)
      call timesec(time_end)
      time_isdf = time_end - time_start

      ! Write results to file
      call h5%initialize_group(h5group, h5group_isdf_wscr_u)
      group = join_paths(h5group, h5group_isdf_wscr_u)
      call h5%write(group, h5ds_indices, indices, [1], [n_isdf])
      call h5%write(group, h5ds_zeta, zeta, [1, 1], [r_grid%number_of_points(), n_isdf])
      call h5%write(group, h5ds_points, r_grid%coordinate_array(), [1, 1], [3, r_grid%number_of_points()])
      call h5%write(group, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1], shape(u_u_isdf))

      call write_info(info_unit, "W_scr unoccupied", time_cvt, time_isdf, has_converged, steps, deviation, n_isdf)

      ! Clean up
      deallocate(indices, zeta, u_u_isdf)
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


  

  !subroutine fastBSE_isdf_qrcp(mpi_env, input, h5file, h5group, info_unit)
  !  type(mpiinfo), intent(inout) :: mpi_env
  !  type(input_type), intent(in) :: input
  !  character(*), intent(in) :: h5file, h5group
  !  integer, intent(in) :: info_unit
!
  !  ! HDF5
  !  type(xhdf5_type) :: h5
  !  character(:), allocatable :: group
!
  !  character(:), allocatable :: seed, bse_type
  !  logical :: calculate_vexc, calculate_wscr
!
  !  integer :: n_r, n_o, n_u, n_k, n_isdf, n_sub, n_sub_input, r_sampling(3), n_combinations, nisdf(3)
  !  integer :: i 
  !  real(dp) :: lattice(3, 3), r_offset(3), epslat, qrcp_eps, tol_qrcp, c_vexc_qrcp, c_wscr_qrcp, time_start, time_end, time_qrcp, time_isdf, r_diag_limit
!
  !  complex(dp), allocatable, target :: u_o(:, :, :), u_u(:, :, :), u_o_isdf(:, :, :), u_u_isdf(:, :, :)
  !  complex(dp), pointer :: u_o_ptr(:, :), u_u_ptr(:, :), u_o_isdf_ptr(:, :), u_u_isdf_ptr(:, :)
  !  
  !  ! QRCP + ISDF
  !  integer, allocatable :: indices(:)
  !  real(dp), allocatable ::  R_diag(:)
  !  complex(dp), allocatable :: M(:, :), tau(:), zeta(:, :)
!
  !  type(regular_grid_type) :: r_grid 
!
!
  !  ! Intialize input data
  !  seed          = input%xs%fastBSE%seed
  !  lattice       = input%structure%crystal%basevect
  !  epslat        = input%structure%epslat
  !  r_offset      = spread(0._dp, 1, 3)
  !  r_sampling    = input%xs%fastBSE%rsampling
  !  r_grid        = setup_unitcell_grid(r_sampling, r_offset, lattice, epslat)
  !  n_r           = r_grid%number_of_points()  
  !  n_o           = input%xs%BSE%nstlbse(2) - input%xs%BSE%nstlbse(1) + 1
  !  n_u           = input%xs%BSE%nstlbse(4) - input%xs%BSE%nstlbse(3) + 1
  !  n_k           = product(input%xs%ngridk)
  !  tol_qrcp    = input%xs%fastBSE%tol_qrcp
  !  c_vexc_qrcp   = input%xs%fastBSE%c_vexc_qrcp
  !  c_wscr_qrcp   = input%xs%fastBSE%c_wscr_qrcp
  !  n_sub_input   = input%xs%fastBSE%n_subsampling_rows
  !  nisdf = input%xs%fastBSE%nisdf
  !  bse_type      = input%xs%BSE%bsetype
  !  call bse_type_to_bool(bse_type, calculate_vexc, calculate_wscr)
!
  !  call set_seed(seed)
!
  !  ! Read wavefunctions
  !  allocate(u_o(n_r, n_o, n_k))
  !  allocate(u_u(n_r, n_u, n_k))
  !  call h5%initialize(h5file, mpi_env%comm)
  !  call h5%read(h5group, h5ds_wfplot, u_o, [1, 1, 1])
  !  call h5%read(h5group, h5ds_wfplot, u_u, [1, n_o + 1, 1])
!
  !
  !  if(calculate_vexc) then
  !    
  !    call timesec(time_start)
!
  !    ! Setup subsampling matrix M
  !    n_combinations = n_k * n_o * n_u
  !    n_sub = min(floor(c_vexc_qrcp * sqrt(1._dp * n_o * n_u * n_k)), n_combinations)
  !    n_isdf = min(nisdf(1), n_o * n_u)
  !    r_diag_limit = 0._dp
!
  !    call setup_subsampling_matrix_samek(mpi_env, n_sub, u_o, u_u, M)
!
  !    ! ISDF grid via QRCP on M
  !    allocate(indices(n_r), source=0)
  !    allocate(tau(n_sub))
  !    call xgeqp3(M, indices, tau)
!
  !    if (tol_qrcp > 0._dp) then
  !      R_diag = abs([(M(i, i), i=1, n_sub)])
  !      r_diag_limit = R_diag(1) * tol_qrcp
  !      n_isdf = count(R_diag >= r_diag_limit)
  !    end if
!
  !    indices = indices(1:n_isdf)
!
  !    call timesec(time_end)
  !    time_qrcp = time_end - time_start
!
  !    call timesec(time_start)
 !
  !    ! Calculate ISDF coefficients zeta
  !    u_o_isdf = u_o(indices, :, :)
  !    u_u_isdf = u_u(indices, :, :)
  !    u_o_ptr(1 : n_r, 1 : n_o * n_k) => u_o
  !    u_u_ptr(1 : n_r, 1 : n_u * n_k) => u_u
  !    u_o_isdf_ptr(1 : n_isdf, 1 : n_o * n_k) => u_o_isdf
  !    u_u_isdf_ptr(1 : n_isdf, 1 : n_u * n_k) => u_u_isdf
  !    call isdf(mpi_env, u_o_ptr, u_u_ptr, u_o_isdf_ptr, u_u_isdf_ptr, indices, zeta)
!
  !    call timesec(time_end)
  !    time_isdf = time_end - time_start
!
  !    ! Write results to file
  !    call h5%initialize_group(h5group, h5group_isdf_vexc)
  !    group = join_paths(h5group, h5group_isdf_vexc)
  !    call h5%write(group, h5ds_indices, indices, [1], [n_isdf])
  !    call h5%write(group, h5ds_zeta, zeta, [1, 1], [r_grid%number_of_points(), n_isdf])
  !    call h5%write(group, h5ds_points, r_grid%coordinate_array(), [1, 1], [3, r_grid%number_of_points()])
  !    call h5%write(group, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1, 1], shape(u_o_isdf))
  !    call h5%write(group, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1, 1], shape(u_u_isdf))
!
  !    call write_info(info_unit, "V_exc", time_qrcp, time_isdf, n_combinations, n_sub, r_diag_limit, n_isdf)
!
  !    ! Clean up
  !    deallocate(indices, M, tau, zeta, u_o_isdf, u_u_isdf)
  !    if(allocated(R_diag)) deallocate(R_diag)
  !    u_o_ptr => null()
  !    u_u_ptr => null()
  !    u_o_isdf_ptr => null()
  !    u_u_isdf_ptr => null()
  !  
  !  end if
!
  !  if (calculate_wscr) then 
!
  !    ! occupied states
!
  !    call timesec(time_start)
!
  !    ! Setup subsampling matrix M
  !    n_combinations = n_k * n_o
  !    n_sub = min(floor(c_wscr_qrcp * sqrt(1._dp * n_combinations)), n_combinations)
  !    n_isdf = min(nisdf(2), n_combinations**2)
  !    r_diag_limit = 0._dp
  !    
  !    call setup_subsampling_matrix_kkp(mpi_env, n_sub, u_o, M)
  !  
  !    ! ISDF grid via QRCP on M
  !    allocate(indices(n_r), source=0)
  !    allocate(tau(n_sub**2))
  !    call xgeqp3(M, indices, tau)
!
  !    if (tol_qrcp > 0._dp) then
  !      R_diag = abs([(M(i, i), i=1, n_sub**2)])
  !      r_diag_limit = R_diag(1) * tol_qrcp
  !      n_isdf = count(R_diag >= r_diag_limit)
  !    end if
!
  !    indices = indices(1:n_isdf)
!
  !    call timesec(time_end)
  !    time_qrcp = time_end - time_start
!
  !    call timesec(time_start)
 !
  !    ! Calculate ISDF coefficients zeta
  !    u_o_isdf = u_o(indices, :, :)
  !    u_o_ptr(1 : n_r, 1 : n_o * n_k) => u_o
  !    u_o_isdf_ptr(1 : n_isdf, 1 : n_o * n_k) => u_o_isdf
  !    call isdf(mpi_env, u_o_ptr, u_o_isdf_ptr, indices, zeta)
!
  !    call timesec(time_end)
  !    time_isdf = time_end - time_start
  !    
  !    ! Write results to file
  !    call h5%initialize_group(h5group, h5group_isdf_wscr_o)
  !    group = join_paths(h5group, h5group_isdf_wscr_o)
  !    call h5%write(group, h5ds_indices, indices, [1], [n_isdf])
  !    call h5%write(group, h5ds_zeta, zeta, [1, 1], [r_grid%number_of_points(), n_isdf])
  !    call h5%write(group, h5ds_points, r_grid%coordinate_array(), [1, 1], [3, r_grid%number_of_points()])
  !    call h5%write(group, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1, 1], shape(u_o_isdf))
!
  !    call write_info(info_unit, "W_scr occupied", time_qrcp, time_isdf, n_combinations**2, n_sub**2, r_diag_limit, n_isdf)
!
  !    ! Clean up
  !    deallocate(indices, M, tau, zeta, u_o_isdf)
  !    if(allocated(R_diag)) deallocate(R_diag)
  !    u_o_ptr => null()
  !    u_o_isdf_ptr => null()
  !    
!
  !    ! unoccupied states
!
  !    call timesec(time_start)
!
  !    ! Setup subsampling matrix M
  !    n_combinations = n_k * n_u
  !    n_sub = min(floor(c_wscr_qrcp * sqrt(1._dp * n_combinations)), n_combinations)
  !    n_isdf = min(nisdf(3), n_combinations**2)
  !    r_diag_limit = 0._dp
!
  !    call setup_subsampling_matrix_kkp(mpi_env, n_sub, u_u, M)
!
  !    ! ISDF grid via QRCP on M
  !    
  !    allocate(indices(n_r), source=0)
  !    allocate(tau(n_sub**2))
  !    call xgeqp3(M, indices, tau)
!
  !    if (tol_qrcp > 0._dp) then
  !      R_diag = abs([(M(i, i), i=1, n_sub**2)])
  !      r_diag_limit = R_diag(1) * tol_qrcp
  !      n_isdf = count(R_diag >= r_diag_limit)
  !    end if
!
  !    indices = indices(1:n_isdf)
!
  !    !call qrcp(mpi_env, n_isdf, n_sub, u_u, indices)
!
  !    call timesec(time_end)
  !    time_qrcp = time_end - time_start
!
  !    call timesec(time_start)
 !
  !    ! Calculate ISDF coefficients zeta
  !    u_u_isdf = u_u(indices, :, :)
  !    u_u_ptr(1 : n_r, 1 : n_u * n_k) => u_u
  !    u_u_isdf_ptr(1 : n_isdf, 1 : n_u * n_k) => u_u_isdf
  !    call isdf(mpi_env, u_u_ptr, u_u_isdf_ptr, indices, zeta)
!
  !    call timesec(time_end)
  !    time_isdf = time_end - time_start
!
  !    ! Write results to file
  !    call h5%initialize_group(h5group, h5group_isdf_wscr_u)
  !    group = join_paths(h5group, h5group_isdf_wscr_u)
  !    call h5%write(group, h5ds_indices, indices, [1], [n_isdf])
  !    call h5%write(group, h5ds_zeta, zeta, [1, 1], [r_grid%number_of_points(), n_isdf])
  !    call h5%write(group, h5ds_points, r_grid%coordinate_array(), [1, 1], [3, r_grid%number_of_points()])
  !    call h5%write(group, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1, 1], shape(u_u_isdf))
!
  !    call write_info(info_unit, "W_scr unoccupied", time_qrcp, time_isdf, n_combinations**2, n_sub**2, r_diag_limit, n_isdf)
!
  !    ! Clean up
  !    deallocate(indices, M, tau, zeta, u_u_isdf)
  !    if(allocated(R_diag)) deallocate(R_diag)
  !    u_u_ptr => null()
  !    u_u_isdf_ptr => null()
  !  end if
!
!  call h5%finalize
!
!  deallocate(u_o, u_u) 
!  
!  contains
!
!  subroutine write_info(info_unit, taskname, time_qrcp, time_isdf, n_combinations, n_sub, tol_times_R_ii, n_isdf)
!    integer, intent(in) :: info_unit
!    character(*), intent(in) :: taskname
!    real(dp) :: time_qrcp, time_isdf
!    integer :: n_combinations 
!    integer :: n_sub
!    real(dp) :: tol_times_R_ii
!    integer :: n_isdf
!    
!    write(info_unit, *)
!    write(info_unit, "(A, A)")     "ISDF + QRCP done for ", taskname
!    write(info_unit, *)
!
!    write(info_unit, "(A, I8)")    "number of combinations:     ", n_combinations
!    write(info_unit, "(A, I8)")    "number of subsampling rows: ", n_sub
!    write(info_unit, *)
!    
!    write(info_unit, "(A, ES21.14)") "lower limit for R_ii: ",  tol_times_R_ii
!    write(info_unit, "(A, I8)")    "number of interpolation points: ",  n_isdf
!    write(info_unit, *)
!
!    write(info_unit, "(A, F15.6)") "Time(QRCP) (s): ", time_qrcp
!    write(info_unit, "(A, F15.6)") "Time(ISDF) (s): ", time_isdf
!    write(info_unit, *)
!    
!  end subroutine write_info
!
!end subroutine fastBSE_isdf_qrcp

  !> Load the periodic part of the wavefunction \(u\) from the hdf5 file.
  subroutine load_u(mpi_env, input, h5file, h5group, u_u, u_o)
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container
    type(input_type), intent(in) :: input
    !> Name of the HDF5 file to read \(u\) from
    character(*), intent(in) :: h5file
    !> Name of the group in the HDF5 file to read \(u\) from
    character(*), intent(in) :: h5group

    complex(dp), allocatable, intent(out) :: u_u(:, :), u_o(:, :)
    
    integer, allocatable :: index_map_u(:, :), index_map_o(:, :)

    integer :: r_sampling(3), n_r, n_k, n_bands, ik, first, last, first_band, last_band
    integer, allocatable :: k_list(:), uo_limits(:, :), n_u(:), n_o(:), shape_u(:)
    complex(dp), allocatable :: u_chunk(:, :, :)

    type(xhdf5_type) :: h5
    character(:), allocatable :: group

    call h5%initialize(h5file, mpi_env%comm)
    group = join_paths(h5group, h5group_wfplot)
    
    call h5%dataset_shape(group, h5ds_u, shape_u, .true.)
    n_r     = shape_u(1)
    n_k     = shape_u(3)
    
    group = join_paths(h5group, h5group_transitions)
    allocate(uo_limits(4, n_k))
    call h5%read(group, h5ds_uo_limits, uo_limits, [1, 1])
    
    n_u = uo_limits(2, :) - uo_limits(1, :) + 1
    n_o = uo_limits(4, :) - uo_limits(3, :) + 1
    
    first_band = uo_limits(3, 1) ! I assume that the lowest occupied and hightes unoccupied band index is constant with ik
    last_band = uo_limits(2, 1)
    n_bands = last_band - first_band + 1

    group = join_paths(h5group, h5group_wfplot)

    if (allocated(u_u)) deallocate(u_u); allocate(u_u(n_r, sum(n_u)))
    if (allocated(u_o)) deallocate(u_o); allocate(u_o(n_r, sum(n_o)))

    allocate(u_chunk(n_r, n_bands, 1))
    do ik=1, n_k
      call h5%read(group, h5ds_u, u_chunk, [1, 1, ik])
      
      first = first_element(n_u, ik)
      last = last_element(n_u, ik)
      u_u(:, first : last) = u_chunk(:, uo_limits(1, ik) : uo_limits(2, ik), 1)

      first = first_element(n_o, ik)
      last = last_element(n_o, ik)
      u_o(:, first : last) = u_chunk(:, uo_limits(3, ik): uo_limits(4, ik), 1)
    end do

    deallocate(u_chunk, uo_limits, n_u, n_o, shape_u)
    
  end subroutine load_u
  

  

end module fastBSE_isdf