!> This module carries global EPH variables that are used across multiple modules.
module eph_variables
  use precision, only: dp
  use mod_kpointset, only: k_set, Gk_set
  use bz_path, only: bz_path_type

  implicit none
  private

  public :: eph_var_init, eph_var_free, eph_var_gen_frequency_grid, eph_var_init_bz_int

  !================================================================================ 
  ! ELECTRON QUANTITIES
  !
  !> \({\bf k}\)-grid on which the electrons and Wannier functions have been calculated
  type(k_set), public :: eph_kset_el
  !> set of \({\bf G+k}\)-vectors for electron wavefunctions
  type(Gk_set), public :: eph_Gkset_el
  !> use of Wannier functions
  logical, public :: eph_use_wannier = .true.
  !> index of first and last electronic state used for Wannierization
  integer, public :: eph_fst, eph_lst
  !> total number of electronic states used for Wannierization
  integer, public :: eph_nst
  !> total number of Wannier functions
  integer, public :: eph_nwf_tot
  !> index of first and last Wannier function to be included in e-ph calculation
  integer, public :: eph_fwf, eph_lwf
  !> number of Wannier functions included in e-ph calculation
  integer, public :: eph_nwf
  !> Fermi energy
  real(dp), public :: eph_efermi
  !> Scissor shift for (un)occupied states
  real(dp), public :: eph_scissor(2)
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! PHONON QUANTITIES
  !
  !> \({\bf q}\)-grid on which the phonons have been calculated
  type(k_set), public :: eph_qset_ph
  !> total number of phonon modes
  integer, public :: eph_nmode_tot
  !> index of first and last phonon mode to be included in e-ph calculation
  integer, public :: eph_fmode, eph_lmode
  !> number of phonon modes included in e-ph calculation
  integer, public :: eph_nmode
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! POLAR MATERIALS
  !
  !> assume polar material
  logical, public :: eph_polar
  !> high-frequency dielectric tensor (clamped-nuclei approximation) \({\bf \epsilon}^\infty\)
  real(dp), public :: eph_dielten(3, 3)
  !> Born-effective charge tensors \({\bf Z}^\ast_\alpha\)
  real(dp), allocatable, public :: eph_borncharge(:,:,:)
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! OTHERS
  !
  !> set of reciprocal space points \({\bf p}\) for which target quantities should be computed 
  !> (can be a grid or a path)
  type(k_set), public :: eph_pset
  type(bz_path_type), public :: eph_pset_path
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! NUMERICAL PARAMETERS
  !
  !> threshold under which phonon frequencies are taken to be zero
  real(dp), parameter, public :: eph_ph_energy_zero = 1e-5_dp
  !> smallest possible finite temperature
  real(dp), parameter, public :: eph_temperature_zero = 1e-16_dp
  !-------------------------------------------------------------------------------- 

contains

  !> Initialize global variables that remain constant during entire EPH calculation.
  subroutine eph_var_init
    use modmpi, only: terminate_if_false
    use phonons_io_util, only: ph_io_read_dielten, ph_io_read_borncharge
    use mod_kpointset, only: generate_k_vectors
    use modinput

    integer :: ip
    logical :: success

    ! check for <wannier> and <phonons> element
    eph_use_wannier = associated(input%properties) 
    if (eph_use_wannier) eph_use_wannier = associated(input%properties%wannier)
    call terminate_if_false( associated(input%phonons), '(eph_var_init) &
      <phonons> element must be present in input file.' )

    ! initialize global variables
    call init0
    call init1

    ! initialize phonons
    call eph_var_init_phonons

    ! initialize electrons
    ! (must happen after phonon initialization in order to keep correct globals for Wannier)
    call eph_var_init_electrons

    ! initialize polar material calculation
    eph_polar = input%phonons%polar
    if (eph_polar) then
      call ph_io_read_dielten( eph_dielten, 'EPSINF.OUT', success )
      call terminate_if_false( success, '(eph_var_init) &
        Failed to read dielectric tensor from file EPSINF.OUT.' )
      call ph_io_read_borncharge( eph_borncharge, 'ZSTAR.OUT', success )
      call terminate_if_false( success, '(eph_var_init) &
        Failed to read Born effective charges from file ZSTAR.OUT.' )
    end if

    ! initialize target point set
    if (associated(input%eph%target)) then
      if (associated(input%eph%target%plot1d)) then
        if (eph_polar) then
          ! slightly shift Gamma point in polar materials in respective direction of path
          eph_pset_path = bz_path_type( eph_kset_el%bvec, input%eph%target%plot1d%path, gamma_offset=1e-6_dp )
        else
          eph_pset_path = bz_path_type( eph_kset_el%bvec, input%eph%target%plot1d%path )
        end if
        call generate_k_vectors( eph_pset, eph_kset_el%bvec, [1, 1, eph_pset_path%num_points], [0.0_dp, 0.0_dp, 0.0_dp], .false., .false. )
        do ip = 1, eph_pset_path%num_points
          eph_pset%vkl(:, ip) = eph_pset_path%points(ip)%coord_lat
          eph_pset%vkc(:, ip) = eph_pset_path%points(ip)%coord_cart
        end do
      else
        call generate_k_vectors( eph_pset, eph_kset_el%bvec, input%eph%target%ngridp, input%eph%target%vploff, input%eph%target%reducep, .false. )
      end if
    else
      call generate_k_vectors( eph_pset, eph_kset_el%bvec, [1, 1, 1], [0.0_dp, 0.0_dp, 0.0_dp], .true., .false. )
    end if
  end subroutine eph_var_init

  !> Free memory from global variables.
  subroutine eph_var_free
    use dfpt_variables, only: dfpt_var_free
    use phonons_variables, only: ph_var_free
    use mod_kpointset, only: delete_k_vectors, delete_Gk_vectors
    use mod_wannier_variables, only: wannier_destroy

    if (allocated(eph_borncharge)) deallocate( eph_borncharge )
    call delete_k_vectors( eph_kset_el )
    call delete_k_vectors( eph_qset_ph )
    call delete_Gk_vectors( eph_Gkset_el )
    call dfpt_var_free
    call ph_var_free
    call wannier_destroy
  end subroutine eph_var_free

  !> Initialize variables of electronic part
  subroutine eph_var_init_electrons
    use mod_wannier, only: wannier_init
    use mod_wannier_variables, only: wf_kset, wf_fst, wf_lst, wf_nwf, wf_efermi
    use modmpi, only: terminate_if_false
    use mod_lattice, only: bvec
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_kpointset, only: generate_k_vectors
    use modinput

    integer :: stype

    ! initialize general variables from Wannier modules
    ! and read Wannier functions from file
    if (eph_use_wannier) then
      input%properties%wannier%do = 'fromfile'
      call wannier_init
    end if

    ! set electron k-grid
    if (eph_use_wannier) then
      eph_kset_el = wf_kset
    else
      stype = input%groundstate%stypenumber
      input%groundstate%stypenumber = 1
      call generate_k_vectors( eph_kset_el, bvec, &
        input%groundstate%ngridk, &
        input%groundstate%vkloff, &
        .false., &
        uselibzint=.false. )
      input%groundstate%stypenumber = stype
    end if

    ! set range of electronic states
    if (eph_use_wannier) then
      where (input%eph%wfrange <= 0) input%eph%wfrange = wf_nwf
      eph_fst = wf_fst
      eph_lst = wf_lst
      eph_nst = eph_lst - eph_fst + 1
    else
      where (input%eph%wfrange <= 0) input%eph%wfrange = nstfv
      eph_fst = max( 1, input%eph%wfrange(1) )
      call terminate_if_false( eph_fst <= nstfv, '(eph_var_init_electrons) &
        Lower bound of `wfrange` must not be greater than total number of states.' )
      eph_lst = min( nstfv, input%eph%wfrange(2) )
      call terminate_if_false( eph_fst <= eph_lst, '(eph_var_init_electrons) &
        Lower bound of `wfrange` must not be greater than upper bound.' )
      eph_nst = eph_lst - eph_fst + 1
    end if

    ! set range of Wannier functions
    if (eph_use_wannier) then
      eph_fwf = max( 1, input%eph%wfrange(1) )
      call terminate_if_false( eph_fwf <= wf_nwf, '(eph_var_init_electrons) &
        Lower bound of `wfrange` must not be greater than total number of Wannier functions.' )
      eph_lwf = min( wf_nwf, input%eph%wfrange(2) )
      call terminate_if_false( eph_fwf <= eph_lwf, '(eph_var_init_electrons) &
        Lower bound of `wfrange` must not be greater than upper bound.' )
      eph_nwf = eph_lwf - eph_fwf + 1
      eph_nwf_tot = wf_nwf
    else
      eph_fwf = 0
      eph_lwf = 0
      eph_nwf = 0
      eph_nwf_tot = 0
    end if

    ! set Fermi energy (all electron energies will be given relative to the Fermi level)
    if (input%eph%efermi /= 0.0_dp) then
      eph_efermi = input%eph%efermi
    else
      eph_efermi = wf_efermi
    end if

    ! set scissor shift
    select case (input%eph%scissor_direction)
      case ('occupied')
        eph_scissor = [input%eph%scissor, 0.0_dp]
      case ('unoccupied')
        eph_scissor = [0.0_dp, input%eph%scissor]
      case default
        eph_scissor = 0.5_dp * input%eph%scissor
    end select
  end subroutine eph_var_init_electrons

  !> Initialize variables of phonon part.
  subroutine eph_var_init_phonons
    use dfpt_variables, only: dfpt_var_init
    use phonons_variables, only: ph_var_init, ph_qset
    use mod_atoms, only: natmtot

    ! initialize general variables of DFPT module
    call dfpt_var_init
    call ph_var_init

    ! set phonon q-grid
    eph_qset_ph = ph_qset

    ! set range and number of phonon modes
    eph_nmode_tot = 3 * natmtot
    eph_fmode = 1
    eph_lmode = eph_nmode_tot
    eph_nmode = eph_nmode_tot
  end subroutine eph_var_init_phonons

  !> Generate grid of frequency points.
  !>
  !> If `energies` is given and `fgrid%type == 'density'`, the grid will be densely sampled around the given energies.
  !> Otherwise, a uniformly sampled grid will be returned.
  function eph_var_gen_frequency_grid( fgrid, energies ) result( freqs )
    use grid_utils, only: linspace, spacing_from_density
    use distributions, only: lorentzian
    use modmpi, only: terminate_if_false
    use modinput, only: freq_grid_type
    !> frequency grid object
    type(freq_grid_type), intent(in) :: fgrid
    !> energies around frequencies should be sampled densely
    real(dp), optional, intent(in) :: energies(:)
    !> frequency grid
    real(dp), allocatable :: freqs(:)
  
    integer, parameter :: nx = 4000 ! number of points for point density
    real(dp), parameter :: ratio = 2.0_dp / (1.0_dp + sqrt( 5.0_dp )) ! mixing ratio between density based and uniform sampling

    integer :: ie

    real(dp), allocatable :: x(:), y(:)

    select case (fgrid%type)
      ! uniform sampling
      case ('uniform')
        freqs = linspace( fgrid%range(1)-fgrid%padding, fgrid%range(2)+fgrid%padding, fgrid%numpoints )
      ! density based sampling
      case ('density')
        call terminate_if_false( present(energies), '(eph_var_gen_frequency_grid): &
          If the grid type is `density`, then `energies` must be provided.' )
        ! generate point density
        x = linspace( fgrid%range(1)-fgrid%padding, fgrid%range(2)+fgrid%padding, nx )
        allocate( y(nx), source=0.0_dp )
        do ie = 1, size(energies)
          y = max( y, lorentzian( x, fgrid%lorentzwidth, energies(ie) ) )
        end do
        y = ratio * y + (1.0_dp - ratio) * sum( y ) / nx
        ! generate frequency grid from point density
        freqs = spacing_from_density( x, y, fgrid%numpoints )
        deallocate( x, y )
    end select
  end function eph_var_gen_frequency_grid

  !> Initialize Brillouin zone integration using tetrahedron integration.
  subroutine eph_var_init_bz_int( ngrid, vloff, pset, tset )
    use mod_kpointset, only: k_set, generate_k_vectors
    use mod_opt_tetra, only: t_set, opt_tetra_init
    use mod_lattice, only: bvec
    use modinput
    !> integration grid
    integer, intent(in) :: ngrid(3)
    !> grid offset in lattice coordinates
    real(dp), intent(in) :: vloff(3)
    !> set of BZ integration points
    type(k_set), intent(out) :: pset
    !> corresponding set of tetrahedra
    type(t_set), intent(out) :: tset

    integer :: stype 

    stype = input%groundstate%stypenumber
    input%groundstate%stypenumber = 1 ! switch off libbzint
    call generate_k_vectors( pset, bvec, ngrid, vloff, .false., uselibzint=.false. )
    call opt_tetra_init( tset, pset, 2, reduce=.false. )
    input%groundstate%stypenumber = stype
  end subroutine eph_var_init_bz_int

end module eph_variables
