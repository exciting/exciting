module fastBSE_write_wfplot
  use precision, only: dp
  use modinput, only: input_type
  use wfplot_nice, only: calculate_wfplot_k_chunk, setup_wfplot_gloabls
  use mod_rgrid, only: rgrid, gen_3d
  use modmpi, only: mpiinfo, distribute_loop
  use grid_utils, only: mesh_1d, concatenate
  use bravais_lattice, only: simple_cubic
  use xhdf5, only: xhdf5_type, abort_if_not_hdf5
  use xfftw, only: abort_if_not_fftw3
  use os_utils, only: join_paths

  implicit none

  private
  public :: fastBSE_write_u

  !>
  character(*), parameter, public :: h5group_wfplot = 'fastBSE_wfplot'
  character(*), parameter, public :: h5ds_u = 'u'
  character(*), parameter, public :: h5ds_r_sampling = 'r_sampling'
  character(*), parameter, public :: h5ds_k_list = 'k_list'
  character(*), parameter, public :: h5ds_band_list = 'h5ds_band_list'
  

  contains

  !> Setup the periodic part of the real space wave functions and write them to an hdf5 file for fast BSE calculations.
  subroutine fastBSE_write_u(h5file, h5group, input, mpi_env)
    use modbse, only: setranges_modxs, select_transitions, koulims, nk_bse
    !> Name of the HDF5 file to write u to.
    character(*), intent(in) :: h5file
    !> Path to the group in to write u to.
    character(*), intent(in) :: h5group
    !> Input file container.
    type(input_type), intent(in) :: input
    !> MPI environment.
    type(mpiinfo), intent(inout) :: mpi_env

    type(rgrid) :: r_grid
    integer, allocatable :: band_list(:), k_list(:)
    integer :: r_sampling(3), offset_u(3), full_shape(3), n_bands,  first, last, nk_local, first_band, last_band
    real(dp) :: box(4, 3)
    complex(dp), allocatable :: u(:, :, :)

    !> Flag for [[gen_3d]] to generate a periodic grid without images.
    integer, parameter :: create_periodic_grid = 1
    integer, parameter :: iqmt = 1

    type(xhdf5_type) :: h5
    character(:), allocatable :: group

    call abort_if_not_fftw3(mpi_env, "Error(fastBSE_write_u): exciting needs to be linked to FFTW3 for running fastBSE.")
    call abort_if_not_hdf5(mpi_env, "Error(fastBSE_write_u): exciting needs to be compiled with HDF5 to run fastBSE module.")

    call setup_wfplot_gloabls(xs_calculation=.true.)
    call readfermi
    call setranges_modxs(iqmt)
    call select_transitions(iqmt, serial=.false.)

    call distribute_loop(mpi_env, nk_bse, first, last)
    k_list = mesh_1d(first, last)
    nk_local = size(k_list)

    first_band = koulims(3, 1) 
    last_band = koulims(2, 1)
    band_list = mesh_1d(first_band, last_band)
    n_bands = size(band_list)

    r_sampling = input%xs%fastBSE%rsampling
    box(1, :) = [0._dp, 0._dp, 0._dp]
    box(2:, :) = simple_cubic(1._dp)
    
    r_grid = gen_3d(r_sampling, box, create_periodic_grid)

    u = calculate_wfplot_k_chunk(r_grid, band_list, k_list, xs_calculation = .true., dephase = .true.)

    offset_u = [1, 1, first]
    full_shape = [r_grid%npt, n_bands, nk_bse]
    
    call h5%initialize(h5file, mpi_env%comm)
    call h5%initialize_group(h5group, h5group_wfplot)
    group = join_paths(h5group, h5group_wfplot)

    call h5%write(group, h5ds_u, u, offset_u, full_shape)
    call h5%write(group, h5ds_r_sampling, r_sampling, [1], [3])
    call h5%write(group, h5ds_k_list, k_list, [first], [nk_bse])
    call h5%write(group, h5ds_band_list, band_list, [1], [n_bands])
    call h5%finalize()
  end subroutine fastBSE_write_u

end module  