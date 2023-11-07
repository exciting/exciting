module fastBSE_transitions
  use precision, only: dp
  use constants, only: zzero
  use math_utils, only: mod1
  use grid_utils, only: first_element, last_element
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

  private
  public :: fastBSE_setup_transitions

  character(*), parameter, public :: h5group_transitions = "fastBSE_transitions"
  character(*), parameter, public :: h5ds_energies = "energies"
  character(*), parameter, public :: h5ds_dmat = "dmat"
  character(*), parameter, public :: h5ds_band_idx = "band_idx"
  character(*), parameter, public :: h5ds_uo_limits = "uo_limits"
  character(*), parameter, public :: h5ds_mask = "mask"


  integer, parameter :: TRUE = 1, FALSE = 0


  contains 

  subroutine fastBSE_setup_transitions(mpi_env, input, h5file, h5group, info_unit)
    use modbse, only: setranges_modxs, select_transitions
    use modbse, only: de, kousize, hamsize, nu_bse_max, no_bse_max, nk_bse, eval0, ofac, smap_rel, koulims, ensortidx
    use m_setup_dmat, only: setup_dmat
    use modmpi, only: mpiglobal
    use constants, only: zi
    use bse_transitions, only: of, ol, uf, ul, tf, tl

    type(mpiinfo), intent(inout) :: mpi_env
    type(input_type), intent(in) :: input
    character(*), intent(in) :: h5file, h5group
    integer, intent(in) :: info_unit

    integer, parameter :: iqmt = 1
    character(*), parameter :: thisname = "fastBSE_setup_transitions"

    real(dp) :: ts1, ts0
    character(:), allocatable :: group 
    type(mpiinfo) :: mpiglobal_save

    integer :: ik, i_transition, i_transition_full, iu, io, nu, no, first_band
    integer, allocatable :: n_o(:), n_u(:), n_uo(:), band_idx(:, :), transition_map_full(:, :), uo_limits(:, :), transition_mask(:)
    complex(dp), allocatable :: dmat(:, :)

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

    ! Write data
    call h5%initialize(h5file, mpi_env%comm)
    call h5%initialize_group(h5group, h5group_transitions)
    group = join_paths(h5group, h5group_transitions)
    call h5%write(group, h5ds_energies, de, [1], shape(de))
    call h5%write(group, h5ds_dmat, dmat, [1, 1], shape(dmat))
    call h5%write(group, h5ds_band_idx, band_idx, [1, 1], shape(band_idx))
    call h5%write(group, h5ds_uo_limits, uo_limits, [1, 1], shape(uo_limits))
    call h5%write(group, h5ds_mask, transition_mask, [1], shape(transition_mask))
    call h5%finalize()

    deallocate(de, dmat)
    if(associated(input%gw)) deallocate(eval0)

  end subroutine fastBSE_setup_transitions

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

end module fastBSE_transitions