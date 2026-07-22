module fastBSE_groundstate_properties
  use precision, only: dp
  use constants, only: zzero
  use math_utils, only: mod1
  use grid_utils, only: first_element, last_element, mesh_1d
  use modmpi, only: mpiinfo, terminate_if_false, distribute_loop
  use modinput, only: input_type

  use xhdf5, only: xhdf5_type, abort_if_not_hdf5
  use xs_hdf5, only: h5ds_wfplot
  use os_utils, only: join_paths
  use xgrid, only: regular_grid_type, setup_unitcell_grid
  use xlapack, only: xgeqp3, qr_column_pivot
  use bse_utils, only: bse_type_to_bool
  use seed_generation, only: set_seed
  use xfftw, only: abort_if_not_fftw3
  use mod_rgrid, only: rgrid, gen_3d
  use m_genfilname, only: genfilname
  use wfplot_nice, only: calculate_wfplot_k_chunk
  use bravais_lattice, only: simple_cubic
  use fastBSE_file_strings

  private
  public :: fastBSE_setup_groundstate_properties, read_wavefunction_u_hdf5, read_transitions_hdf5




  integer, parameter :: TRUE = 1, FALSE = 0


  contains 


  !> Read all relevant data according to the transitions from an HDF5 file.
  subroutine read_transitions_hdf5(mpi_env, h5file, h5path, transition_energies, matrix_elements, band_index, transition_mask)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Name of the HDF5 file
    type(xhdf5_type), intent(inout) :: h5file
    !> Name of the group in the HDF5 file
    character(*), intent(in) :: h5path
    !> Independent particle transition energies \(\epsilon_{u \mathbf k} - \epsilon_{o \mathbf k)\)
    real(dp), allocatable, intent(out) :: transition_energies(:)
    !> Renormalized matrix elements 
    !> \( \frac {\langle o \mathbf k | \mathbf p | u \mathbf k \rangle }{epsilon_{u \mathbf k} - \epsilon_{o \mathbf k)}\)
    complex(dp), allocatable, intent(out) :: matrix_elements(:, :)
    !> Band index lookup table
    integer, allocatable, intent(out) :: band_index(:, :)
    !> Integer mask for allowed transitions
    integer, allocatable, intent(out) :: transition_mask(:)

    integer, allocatable :: shape_transition_energies(:), shape_matrix_elements(:), &
                            shape_band_index(:), shape_transition_mask(:)
    integer :: n_transitions, n_k, ngridk(3)
    character(:), allocatable :: group 

    group = join_paths(h5path, groundstate_properties_group)

    call h5file%dataset_shape(group, transition_energies_dataset, shape_transition_energies)
    call h5file%dataset_shape(group, matrix_elements_dataset, shape_matrix_elements, complex_dataset = .true.)
    call h5file%dataset_shape(group, band_index_dataset, shape_band_index)
    call h5file%dataset_shape(group, transition_mask_dataset, shape_transition_mask)
    call h5file%read(group, ngridk_dataset, ngridk, [1])

    n_transitions = shape_transition_energies(1)
    n_k = product(ngridk)

    call terminate_if_false(mpi_env, shape_matrix_elements(1) == n_transitions, &
            'read_transitions_hdf5: For dataset ' // matrix_elements_dataset // ' dimension 1 is not &
            the same as the size of dataset ' // transition_energies_dataset // ' (n_transitions) as expected.')

    call terminate_if_false(mpi_env, shape_matrix_elements(2) == 3, &
            'read_transitions_hdf5: For dataset ' // matrix_elements_dataset // ' dimension 2 is not 3.')

    call terminate_if_false(mpi_env, shape_band_index(1) ==  6, &
            'read_transitions_hdf5: For dataset ' // band_index_dataset // ' dimension 1 is not 6.')

    call terminate_if_false(mpi_env, shape_band_index(2) == n_k, &
            'read_transitions_hdf5: For dataset ' // band_index_dataset // ' dimension 2 is not the same &
            as the number of k-points (n_k).')

    call terminate_if_false(mpi_env, shape_transition_mask(1) == n_transitions, &
            'read_transitions_hdf5: For dataset ' // transition_mask_dataset // ' dimension 1 is not &
            the same as the size of dataset ' // transition_energies_dataset // ' (n_transitions) as expected.')

    if(allocated(transition_energies)) deallocate(transition_energies)
    allocate(transition_energies(n_transitions))
    call h5file%read(group, transition_energies_dataset, transition_energies)

    if(allocated(matrix_elements)) deallocate(matrix_elements)
    allocate(matrix_elements(n_transitions, 3))
    call h5file%read(group, matrix_elements_dataset, matrix_elements)

    if(allocated(band_index)) deallocate(band_index)
    allocate(band_index(6, n_k))
    call h5file%read(group, band_index_dataset, band_index)

    if(allocated(transition_mask)) deallocate(transition_mask)
    allocate(transition_mask(n_transitions))
    call h5file%read(group, transition_mask_dataset, transition_mask)
  end subroutine 

  !> Load the periodic part of the wavefunction \(u\) from the hdf5 file.
  subroutine read_wavefunction_u_hdf5(mpi_env, input, h5file, h5group, u_u, u_o)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Input file container
    type(input_type), intent(in) :: input
    !> Name of the HDF5 file to read \(u\) from
    type(xhdf5_type), intent(inout) :: h5file
    !> Name of the group in the HDF5 file to read \(u\) from
    character(*), intent(in) :: h5group

    complex(dp), allocatable, intent(out) :: u_u(:, :), u_o(:, :)
    
    integer, allocatable :: index_map_u(:, :), index_map_o(:, :)

    integer :: n_r, n_k, n_bands, ik, first, last, first_band, last_band
    integer, allocatable :: k_list(:), uo_limits(:, :), n_u(:), n_o(:), shape_u(:)
    complex(dp), allocatable :: u_chunk(:, :, :)
    character(:), allocatable :: group

    group = join_paths(h5group, groundstate_properties_group)
    call h5file%dataset_shape(group, u_dataset, shape_u, .true.)
    n_r     = shape_u(1)
    n_k     = shape_u(3)    
    allocate(uo_limits(4, n_k))
    call h5file%read(group, uo_limits_dataset, uo_limits, [1, 1])
    n_u = uo_limits(2, :) - uo_limits(1, :) + 1
    n_o = uo_limits(4, :) - uo_limits(3, :) + 1
    
    first_band = uo_limits(3, 1) ! I assume that the lowest occupied and highest unoccupied band index is constant with ik
    last_band = uo_limits(2, 1)
    n_bands = last_band - first_band + 1

    if (allocated(u_u)) deallocate(u_u)
    allocate(u_u(n_r, sum(n_u)))

    if (allocated(u_o)) deallocate(u_o)
    allocate(u_o(n_r, sum(n_o)))

    allocate(u_chunk(n_r, n_bands, 1))
    do ik=1, n_k
      call h5file%read(group, u_dataset, u_chunk, [1, 1, ik])
      
      first = first_element(n_u, ik)
      last = last_element(n_u, ik)
      u_u(:, first : last) = u_chunk(:, uo_limits(1, ik) : uo_limits(2, ik), 1)

      first = first_element(n_o, ik)
      last = last_element(n_o, ik)
      u_o(:, first : last) = u_chunk(:, uo_limits(3, ik): uo_limits(4, ik), 1)
    end do

    deallocate(u_chunk, uo_limits, n_u, n_o, shape_u) 
  end subroutine read_wavefunction_u_hdf5

  !> Setup all relevant data according to transitions for **fastBSE**.
  subroutine fastBSE_setup_groundstate_properties(mpi_env, input, h5file, h5group, info_unit)
    use modbse, only: setranges_modxs, select_transitions
    use modbse, only: de, kousize, hamsize, nu_bse_max, no_bse_max, nk_bse, eval0, ofac, smap_rel, koulims, ensortidx
    use m_setup_dmat, only: setup_dmat
    use modmpi, only: mpiglobal
    use constants, only: zi
    use bse_transitions, only: of, ol, uf, ul, tf, tl
    use mod_lattice, only: avec
    use modxs, only: evalsv0

    type(mpiinfo), intent(in) :: mpi_env
    type(input_type), intent(in) :: input
    character(*), intent(in) :: h5file, h5group
    integer, intent(in) :: info_unit

    !> Flag for [[gen_3d]] to generate a periodic grid without images.
    integer, parameter :: create_periodic_grid = 1
    !> Index of the momentum transition. Yet only \(\Gamma\)-point is supported.
    integer, parameter :: iqmt = 1
    !> Name of this routine
    character(*), parameter :: thisname = "fastBSE_groundstate_properties"

    real(dp) :: ts1, ts0
    character(:), allocatable :: group 
    type(mpiinfo) :: mpiglobal_save

    integer :: ik, i_transition, i_transition_full, iu, io, nu, no, first_band, last_band
    integer, allocatable :: n_o(:), n_u(:), n_uo(:), band_idx(:, :), transition_map_full(:, :), uo_limits(:, :), transition_mask(:)
    complex(dp), allocatable :: dmat(:, :)
    real(dp), allocatable :: eigen_energies(:, :)

    type(rgrid) :: r_grid
    integer, allocatable :: band_list(:), k_list(:)
    integer :: ngridr(3), offset_u(3), full_shape(3), n_bands,  first, last, nk_local
    real(dp) :: box(4, 3)
    complex(dp), allocatable :: u(:, :, :)

    type(xhdf5_type) :: h5

    call abort_if_not_fftw3(mpi_env, "Error(fastBSE_write_u): exciting needs to be linked to FFTW3 for running fastBSE.")
    call abort_if_not_hdf5(mpi_env, "Error(fastBSE_write_u): exciting needs to be compiled with HDF5 to run fastBSE module.")

    ! Save mpiglobal and set it to mpi_env
    mpiglobal_save = mpiglobal
    mpiglobal = mpi_env 

    ! Setup exciting globals
    call timesec(ts0)
    call init0
    call init1
    call xssave0
    call init2
    call timesec(ts1)
    call readfermi
    call setranges_modxs(iqmt)
    call genfilname(iqmt=iqmt, setfilext=.true.)
    write(info_unit, '("Info(",a,"): Init time: ", f12.6)') trim(thisname), ts1 - ts0

    if(associated(input%gw)) call load_qp_energies(info_unit)

    ! Setup transition energies
    call select_transitions(iqmt, serial=.false.)
    !Setup dipole matrix elements
    allocate(dmat(hamsize, 3))
    call setup_dmat(dmat)

    first_band = koulims(3, 1)
    uo_limits = koulims - first_band + 1 
    n_u = uo_limits(2, :) - uo_limits(1, :) + 1
    n_o = uo_limits(4, :) - uo_limits(3, :) + 1
    n_uo = kousize

    ! Setup look up table for band limits per k-point
    allocate(band_idx(6, nk_bse))
    do ik=1, nk_bse
      band_idx(uf, ik) = first_element(n_u, ik)
      band_idx(of, ik) = first_element(n_o, ik)
      band_idx(tf, ik) = first_element(n_uo, ik)
      band_idx(ul, ik) = last_element(n_u, ik)
      band_idx(ol, ik) = last_element(n_o, ik)
      band_idx(tl, ik) = last_element(n_uo, ik)
    end do

    allocate(transition_map_full(3, sum(n_o * n_u)))
    i_transition = 1
    do ik=1, nk_bse
      nu = n_u(ik)
      no = n_o(ik)
      do io=1, no
        do iu=1, nu
          transition_map_full(:, i_transition) = [iu, io, ik]
          i_transition = i_transition + 1
        end do 
      end do
    end do

    allocate(transition_mask(sum(n_o * n_u)), source=1)
    
    i_transition_full = 1
    do i_transition=1, hamsize
      do while(any(smap_rel(:, i_transition) /= transition_map_full(:, i_transition_full)))
        transition_mask(i_transition_full) = 0
        i_transition_full = i_transition_full + 1
      end do 
      i_transition_full = i_transition_full + 1
    end do

    ! Reset mpiglobal
    mpiglobal = mpiglobal_save


    ! Calculate real
    call distribute_loop(mpiglobal, nk_bse, first, last)
    k_list = mesh_1d(first, last)
    nk_local = size(k_list)

    first_band = koulims(3, 1) 
    last_band = koulims(2, 1)
    band_list = mesh_1d(first_band, last_band)
    n_bands = size(band_list)

    eigen_energies = evalsv0(first_band : last_band, first : last)

    ngridr = input%xs%fastBSE%ngridr
    box(1, :) = [0._dp, 0._dp, 0._dp]
    box(2:, :) = simple_cubic(1._dp)
    
    r_grid = gen_3d(ngridr, box, create_periodic_grid)

    u = calculate_wfplot_k_chunk(r_grid, band_list, k_list, xs_calculation = .true., dephase = .true.)

    offset_u = [1, 1, first]
    full_shape = [r_grid%npt, n_bands, nk_bse]

    ! Write data
    call h5%initialize(h5file, mpiglobal)
    call h5%initialize_group(h5group, groundstate_properties_group)
    group = join_paths(h5group, groundstate_properties_group)
    call h5%write(group, transition_energies_dataset, de)
    call h5%write(group, matrix_elements_dataset, dmat)
    call h5%write(group, band_index_dataset, band_idx)
    call h5%write(group, uo_limits_dataset, uo_limits)
    call h5%write(group, transition_mask_dataset, transition_mask)
    call h5%write(group, u_dataset, u, offset_u, full_shape)
    call h5%write(group, ngridr_dataset, ngridr, [1], [3])
    call h5%write(group, k_list_dataset, k_list, [first], [nk_bse])
    call h5%write(group, band_list_dataset, band_list)
    call h5%write(group, ngridk_dataset, input%xs%ngridk)
    call h5%write(group, ngridr_dataset, ngridr)
    call h5%write(group, lattice_vectors_dataset, avec)
    call h5%write(group, eigen_energies_dataset, eigen_energies, [1, first], [n_bands, nk_bse])
    call h5%finalize()

    deallocate(de, dmat, eigen_energies)
    if(associated(input%gw)) deallocate(eval0)
  end subroutine fastBSE_setup_groundstate_properties


  !> Load quasi particle energies from a GW calculation.
  subroutine load_qp_energies(info_unit)
    ! Globals 
    use modbse, only: eval0
    use mod_eigenvalue_occupancy, only: evalsv, nstsv
    use modxs, only: vkl0, evalsv0
    use mod_kpoint, only: nkptnr
    use mod_symmetry, only: nsymcrys
    use mod_wannier_bse, only: wfbse_usegwwannier, wfbse_init, wfbse_ordereval, wfbse_eval

    !> Info file unit for user output.
    integer, intent(in) :: info_unit
    
    character(*), parameter :: thisname='load_qp_energies'

    integer(4) :: nsymcrys_save

    ! Save KS eigenvalues of the k-grid to use them later for renormalizing PMAT
    if(allocated(eval0)) deallocate(eval0)
    allocate(eval0(nstsv, nkptnr))
    eval0=evalsv0

    ! Read QP Fermi energies and eigenvalues from file
    ! NOTE: QP evals are shifted by -efermi-eferqp with respect to KS evals
    ! NOTE: getevalqp sets mod_symmetry::nsymcrys to 1
    ! NOTE: getevalqp needs the KS eigenvalues as input
    if( wfbse_usegwwannier()) then
      call wfbse_init
      call wfbse_ordereval
      evalsv = wfbse_eval
    else
      nsymcrys_save = nsymcrys
      !call checkevalqp('EVALQP.OUT', nkptnr, vkl0, evalsv)
      call getevalqp('EVALQP.OUT', nkptnr, vkl0, evalsv)
      nsymcrys = nsymcrys_save
    end if

    ! Set k and k'=k grid eigenvalues to QP energies
    evalsv0=evalsv

    write(info_unit,'("Info(",a,"): Quasi particle energies are read from EVALQP.OUT")') thisname
    if( wfbse_usegwwannier()) then
      write(info_unit,'("Info(",a,"): Wannier interpolation was employed.")') thisname
    end if
  end subroutine 

end module fastBSE_groundstate_properties