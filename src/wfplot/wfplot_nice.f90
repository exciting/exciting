module wfplot_nice
  use precision, only: dp
  use asserts, only: assert
  use constants, only: zone
  use modinput, only: input
  use modmpi, only: terminate_if_false
  use mod_rgrid, only: rgrid, calc_zdata_rgrid
  use m_genfilname
  use modxs, only: iqmtgamma
  use grid_utils, only: phase

  private
  public :: setup_wfplot_gloabls, setup_wfplot_k, calculate_wfplot, calculate_wfplot_k_chunk

  contains 

  !> Setup globals and loop independent prerequirements for calculating the wave function on a real
  !> space r_grid.
  !> The routines sets `input%groundstate%lradstep` to `1`!
  subroutine setup_wfplot_gloabls(xs_calculation)
    
    use mod_atoms, only: natmtot
    use mod_Gkvector, only: ngkmax
    use mod_APW_LO, only: apwordmax
    use mod_muffin_tin, only: lmmaxapw
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv, nstsv

   
    !> If the wfplot calculation is in the context of a xs calculation set this flag to true to 
    !> initialize the required global variables (`call init2()`). See [[xs_calculation_default]]
    !> for the default value.
    logical, intent(in) :: xs_calculation



    ! initialise universal variables
    input%groundstate%lradstep = 1
    
    call init0()
    call init1()

    if(xs_calculation) then 
      call init2()
      call xssave0
      call genfilname(iqmt=iqmtgamma, setfilext=.true.)
    else
      ! read the density and potentials from file
      call readstate
      ! find the new linearisation energies
      call linengy
      ! generate the APW radial functions
      call genapwfr
      ! generate the local-orbital radial functions
      call genlofr
    end if

  end subroutine setup_wfplot_gloabls

  !> Setup arrays that are needed for calculating the wfplot and that only depend on the \(\mathbf{k}\)-point but not
  !> on the state number.
  subroutine setup_wfplot_k(ik, evecfv, evecsv, apwalm)
    ! Globals
    
    use mod_kpoint, only: vkl
    use mod_Gkvector, only: ngk, vgkl, gkc, tpgkc, sfacgk

    !> Index of the k point
    integer, intent(in) :: ik
    !> First variational eigenvectors for \(\mathbf{k}\)-point ik
    complex(dp), intent(inout) :: evecfv(:, :)
    !> Second variational eigenvectors for \(\mathbf{k}\)-point ik
    complex(dp), intent(inout) :: evecsv(:, :)
    !> APW matching coefficients for \(\mathbf{k}\)-point ik
    complex(dp), intent(inout) :: apwalm(:, :, :, :)

    ! get the eigenvectors and values from file
    call getevalsv(vkl(:,ik), evecsv)
    call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfv)
    call getevecsv(vkl(:,ik), evecsv)

    ! find the matching coefficients
    call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), sfacgk(:,:,1,ik), apwalm)
  end subroutine setup_wfplot_k

  !> Calculate the wave function on a real space r_grid for a \(\mathbf{k}\)-point and a state,
  !> specified by their indices.
  !>
  !> CAUTION: This routine ist highly dependent on globals. It requires to call [[setup_wfplot_globals]] beforehand.
  function calculate_wfplot(ik, ist, r_grid, evecfv, evecsv, apwalm) result(zdata)

    use mod_atoms, only: natmtot
    use mod_Gkvector, only: ngkmax
    use mod_APW_LO, only: apwordmax
    use mod_muffin_tin, only: lmmaxapw
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv, nstsv

    use mod_muffin_tin, only: lmmaxapw, nrmtmax
    use mod_atoms, only: natmtot
    use mod_spin, only: nspinor
    use mod_eigenvalue_occupancy, only: nstsv
    use mod_kpoint, only: nkpt
    use mod_Gvector, only: ngrtot

    !> Index of the \(\mathbf{k}\)-point for which the wavefunction plot is calculated
    integer, intent(in) :: ik
    !> Index of the state for which the wavefunction plot is calculated
    integer, intent(in) :: ist
    !> Meta data for the r_grid
    type(rgrid), intent(in) :: r_grid
    !> First variational eigenvectors for \(\mathbf{k}\)-point ik
    complex(dp), intent(in) :: evecfv(:, :)
    !> Second variational eigenvectors for \(\mathbf{k}\)-point ik
    complex(dp), intent(in) :: evecsv(:, :)
    !> APW matching coefficients for \(\mathbf{k}\)-point ik
    complex(dp), intent(in) :: apwalm(:, :, :, :)

    complex(dp), allocatable :: zdata(:)

    complex(dp), allocatable ::  wfmt(:,:,:,:,:), wfir(:,:,:)

    call assert(1 <= ik, 'k-point index ist smaller than 1 (ik < 1).')
    call assert(ik <= nkpt, 'k-point index ist larger than the number of k-points (ik > nkpt).')
    call assert(1 <= ist, 'State index ist smaller than 1 (ist < 1).')
    call assert(ist <= nstsv, 'State index ist larger than the number states (ik > nstsv).')

    ! calculate the wavefunctions for all states
    allocate(wfmt(lmmaxapw, nrmtmax, natmtot, nspinor, nstsv))
    allocate(wfir(ngrtot, nspinor, nstsv))
    call genwfsv_new(ik, ist, ist, apwalm, evecfv, evecsv, wfmt, wfir)

    allocate(zdata(r_grid%npt))
    call calc_zdata_rgrid(r_grid, ik, wfmt(:, :, :, 1, ist), wfir(:, 1, ist), zdata)
  end function calculate_wfplot

  !> Calculate the `wfplot` for the interval of \(\mathbf{k}\)-points \([\mathbf{k}_{first}, \mathbf{k}_{last}]\).
  !> [[setup_wfplot_globals]] must be called before this routine.
  function calculate_wfplot_k_chunk(r_grid, band_list, k_list, xs_calculation, dephase) result(wfplot_chunk)
    use mod_kpoint, only: vkl, nkpt
    use mod_atoms, only: natmtot
    use mod_Gkvector, only: ngkmax
    use mod_APW_LO, only: apwordmax
    use mod_muffin_tin, only: lmmaxapw
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    !> \(\mathbf{r}\)-grid to calculate the wave functions on
    type(rgrid), intent(in) :: r_grid
    !> List of band indices to calculate the wave functions for
    integer, intent(in) :: band_list(:)
    !> List of \(\mathbf{k}\)-points to calculate the wave functions for
    integer, intent(in) :: k_list(:)
    !> Flag to set true for xs calculation. Makes sure that consistant file names are used.
    logical, intent(in) :: xs_calculation
    !> If true, the output is multiplied with the inverse Bloch phase to make it periodic.
    logical, intent(in) :: dephase

    complex(dp), allocatable :: wfplot_chunk(:, :, :) ! Data to write N_r x N_st x N_k

    logical :: xs_calculation_local

    integer :: i, ist, ik, N_st, N_k, N_r
    complex(dp), allocatable :: apwalm(:, :, :, :), evecfv(:, :), evecsv(:, :)

    call assert(minval(band_list) > 0, 'minval(band_list) <= 0.')
    call assert(maxval(band_list) <= nstsv, 'maxval(band_list) <= nstsv.')
    call assert(minval(k_list) > 0, 'minval(k_list) <= 0.')
    call assert(maxval(k_list) <= nkpt, 'maxval(k_list) <= nkpt.')

    N_r = r_grid%npt
    N_st = size(band_list)
    N_k = size(k_list)

    allocate(wfplot_chunk(N_r, N_st, N_k))
    allocate(evecfv(nmatmax, nstfv))
    allocate(evecsv(nstsv, nstsv))
    allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))

    do ik = 1, N_k
      call setup_wfplot_k(k_list(ik), evecfv, evecsv, apwalm)
      do ist = 1, N_st
        wfplot_chunk(:, ist, ik) = calculate_wfplot(k_list(ik), band_list(ist), r_grid, evecfv, evecsv, apwalm)
      end do
      if (dephase) then
        wfplot_chunk(:, :, ik) = wfplot_chunk(:, :, ik) * spread(phase(r_grid%vpl, vkl(:, k_list(ik))), 2, N_st)
      end if  
    end do
  end function calculate_wfplot_k_chunk

end module wfplot_nice
