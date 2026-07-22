module dfthalf
  use asserts, only: assert
  use mod_potential_and_density, only: vhalfir, vhalfmt
  use modmpi, only : terminate_if_false
  use precision, only: dp, i32
  use vhalfinit, only: initialize_vhalf
  use to_char_conversion, only: to_char
#include "asserts.fpp"
  implicit none

  private

  public :: allocate_vhalf_global_arrays, &
            obtain_vhalf_potential, &
            run_dft_half_nscf

  type, private :: dft_half_shell_parameters
  private
    integer(i32) :: number = 0
    real(dp) :: ionization = 0._dp
    real(dp) :: ampl = 1._dp
    real(dp) :: cut = 0._dp
    integer(i32) :: exponent = 8
  contains
  private
    procedure, pass :: parse_input => dft_half_shell_parameters_parse_input
    procedure, pass :: sanity_check => dft_half_shell_parameters_sanity_check
  end type

  type, private :: dft_half_species_parameters
  private
    !> Maximum number of shells for the given species
    integer(i32) :: max_shells = 0
    !> Array of shell parameters for the given species
    type(dft_half_shell_parameters), allocatable :: shell_array(:)
  contains
  private
    procedure, pass :: parse_input => dft_half_species_parameters_parse_input
    procedure, pass :: sanity_check => dft_half_species_parameters_sanity_check
  end type

  type, public :: dft_half_parameters
  private
    logical :: printVSfile = .false.
    type(dft_half_species_parameters), allocatable :: species_array(:)
  contains
    procedure, pass :: parse_input => dft_half_parameters_parse_input
    procedure, pass :: sanity_check => dft_half_parameters_sanity_check
  end type

contains

pure subroutine dft_half_shell_parameters_parse_input( this, shell_input, max_shells )
  use modinput, only: shell_type
  class(dft_half_shell_parameters), intent(inout) :: this
  !> Input shell parameters
  type(shell_type), intent(in) :: shell_input
  !> Maximum number of shells (for the a given species)
  integer(i32), intent(in) :: max_shells

  this%number = shell_input%number
  this%ionization = shell_input%ionization
  this%ampl = shell_input%ampl
  this%cut = shell_input%cut
  this%exponent = shell_input%exponent
  ! Fix non-positive shell numbers
  if( this%number <= 0 ) this%number = max_shells + this%number
end subroutine

pure logical function dft_half_shell_parameters_sanity_check( this, max_shells ) result( ok )
  class(dft_half_shell_parameters), intent(in) :: this
  !> Maximum number of shells (for the a given species)
  integer(i32), intent(in) :: max_shells

  ok = ( this%number > 0 ) .and. ( this%number <= max_shells )
end function

pure subroutine dft_half_species_parameters_parse_input( this, shell_array, max_shells )
  use modinput, only: shell_type_array
  class(dft_half_species_parameters), intent(inout) :: this
  !> Array of shell parameters for the given species (from the input file)
  type(shell_type_array), intent(in) :: shell_array(:)
  !> Maximum number of shells (for the a given species)
  integer(i32), intent(in) :: max_shells

  integer(i32) :: i, n_shells

  n_shells = size( shell_array )
  this%max_shells = max_shells
  allocate( this%shell_array( n_shells ) )
  do i = 1, n_shells
    call this%shell_array(i)%parse_input( shell_array(i)%shell, max_shells )
  end do
end subroutine

pure logical function dft_half_species_parameters_sanity_check( this ) result( ok )
  class(dft_half_species_parameters), intent(in) :: this

  integer(i32) :: i

  ok = .true.
  do i = 1, size( this%shell_array )
    ok = ok .and. this%shell_array(i)%sanity_check( this%max_shells )
  end do
end function

pure subroutine dft_half_parameters_parse_input( this, dft_half_input, species_array, max_shells_array )
  use modinput, only: species_type_array, dfthalf_type
  class(dft_half_parameters), intent(inout) :: this
  !> DFT-1/2 input parameters
  type(dfthalf_type), intent(in) :: dft_half_input
  !> Array with the species parameters (from the input file)
  type(species_type_array), pointer, intent(in) :: species_array(:)
  !> Array with the maximum number of shells for each species
  integer(i32), intent(in) :: max_shells_array(:)

  integer(i32) :: i, n_species

  this%printVSfile = dft_half_input%printVSfile
  n_species = size( species_array )
  allocate( this%species_array(n_species) )
  do i = 1, n_species
    if( associated( species_array(i)%species%dfthalfparam ) ) then
      call this%species_array(i)%parse_input( species_array(i)%species%dfthalfparam%shellarray, max_shells_array(i) )
    end if
  end do
end subroutine

subroutine dft_half_parameters_sanity_check( this )
  class(dft_half_parameters), intent(in) :: this

  integer(i32) :: i

  do i = 1, size( this%species_array )
    if( allocated( this%species_array(i)%shell_array ) ) then
      call terminate_if_false( this%species_array(i)%sanity_check( ), &
        'Error concerning the shell%number of species '// to_char( i ) // &
        'It should not be negative or higher than the number of total shells' )
    end if
  end do
end subroutine

!> Allocate global arrays needed for DFT-1/2 calculations
subroutine allocate_vhalf_global_arrays( n_r_IR, n_r_MT, lm_max, n_atoms )
  !> Number of radial points in the interstitial region
  integer(i32), intent(in) :: n_r_IR
  !> Number of radial points in the muffin-tin region
  integer(i32), intent(in) :: n_r_MT
  !> Maximum number of angular momentum components used for the expansion of the potential 
  !> inside muffin-tin spheres
  integer(i32), intent(in) :: lm_max
  !> Total number of atoms in the system
  integer(i32), intent(in) :: n_atoms

  if( allocated( vhalfir ) ) deallocate( vhalfir )
  if( allocated( vhalfmt ) ) deallocate( vhalfmt )
  allocate( vhalfir(n_r_IR) )
  allocate( vhalfmt(lm_max, n_r_MT, n_atoms) )
end subroutine

!> Obtain the DFT-1/2 potential
subroutine obtain_vhalf_potential( xctype, xcgrad, dirac_eq, point_nucleus, params )
  use mod_atoms, only: nspecies, speval, spk, spl, spn, spnr, spnst, spnrmax, spnstmax, spocc, spr, spvr, spzn
  use mod_muffin_tin, only: nrmt
  use modmpi, only : mpiglobal

  !> Type of exchange-correlation functional
  integer(i32), intent(in) :: xctype
  !> Gradient of the exchange functional
  integer(i32), intent(in) :: xcgrad
  !> If `.true.` solve the Dirac equation
  logical, intent(in) :: dirac_eq
  !> `.true.` if the nucleus must be treated as a point particle
  logical, intent(in) :: point_nucleus
  !> Object containing DFT-1/2 parameters
  type(dft_half_parameters), intent(in) :: params

  integer(i32) :: is, i, n_shell, nr, xctypearray(3)
  real(dp), allocatable :: cut_function(:), rwf(:, :, :), new_occupations(:)
  real(dp), allocatable :: rho(:), v_KS(:), delta_v(:), VS_untrimmed(:)
  real(dp), allocatable :: VS_trimmed(:, :)

  xctypearray(1:3) = xctype
  allocate( rwf(spnrmax, 2, spnstmax), rho(spnrmax) )
  allocate( VS_trimmed(spnrmax, nspecies), source = 0._dp ) 
  allocate( v_KS(spnrmax), delta_v(spnrmax), VS_untrimmed(spnrmax), cut_function(spnrmax) ) 
  allocate( new_occupations(spnstmax) )
  do is = 1, nspecies
    if( allocated( params%species_array(is)%shell_array ) ) then 
      n_shell = size( params%species_array(is)%shell_array )
      if( params%printVSfile ) VS_untrimmed = 0._dp
      do i = 1, n_shell
        call set_new_occupations( is, params%species_array(is)%shell_array(i), &
          spnst(is), spocc(:, is), new_occupations )
        ! Calculate the ionized atom
        call atom ( point_nucleus, spzn(is), spnst(is), spn(:, is), spl(:, is), &
          spk(:, is), new_occupations, xctypearray, xcgrad, spnr(is), spr(:, is), &
          speval(:, is), rho, v_KS, rwf, nrmt(is), dirac_eq )
        delta_v = v_KS - spvr(:, is)
        call obtain_cut_function( spr(:, is), params%species_array(is)%shell_array(i), cut_function )
        VS_trimmed(:, is) = VS_trimmed(:, is) + cut_function*delta_v
        if( params%printVSfile ) VS_untrimmed = VS_untrimmed + delta_v
      end do
      if( params%printVSfile .and. mpiglobal%is_root ) then
        nr = spnr(is)
        call write_VS_to_file( is, params%species_array(is)%shell_array, &
          spr(1:nr, is), spvr(1:nr, is), VS_untrimmed(1:nr), VS_trimmed(1:nr, is) )
      end if
    end if
  end do
  ! Now, expand the V_S potential in the MT and in the interstitial part
  call initialize_vhalf( VS_trimmed )
end subroutine

subroutine run_dft_half_nscf( )
  use exciting_mpi, only: xmpi_gatherv
  use modmpi, only: distribute_loop, mpiglobal
  use mod_atoms, only: natmtot
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use mod_gen_lo, only: genlofr
  use mod_getoccsv, only: getoccsv
  use mod_gkvector, only: gkmax
  use mod_gvector, only: intgv
  use mod_kpointset, only: generate_G_vectors, generate_Gk_vectors, generate_k_vectors, &
    G_set, Gk_set, k_set
  use mod_lattice, only: bvec
  use modinput, only: input
  
  integer(i32) :: ik, first_kpt, last_kpt, n_states, n_states_sv, l_max
  real(dp), allocatable :: ks_eigenvalues(:, :), occupations(:, :)
  real(dp), allocatable :: buffer_s_alpha_matrix(:, :), s_alpha_matrix(:, :)
  real(dp), allocatable, target :: buffer_band_character(:, :, :, :), buffer_aux(:, :)
  real(dp), contiguous, pointer :: buffer_pointer(:, :), band_character(:, :, :, :)
  type(k_set) :: kset
  type(G_set) :: Gset
  type(Gk_set) :: Gkset
  logical :: print_band_character

  call init0
  call init1
  call init2

  call readstate()        ! read the density and potentials from file
  call gencore()          ! generate the core wavefunctions and densities
  call genmeffig()
  call linengy()          ! find the new linearization energies
  call genapwfr()         ! generate the APW radial functions
  call genlofr()          ! generate the local-orbital radial functions
  call olprad()           ! compute the overlap radial integrals

  call generate_k_vectors( kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, &
    input%groundstate%reducek, .false. )
  call distribute_loop( mpiglobal, kset%nkpt, first_kpt, last_kpt )
  n_states = nstfv
  n_states_sv = nstsv

  call generate_G_vectors( Gset, bvec, intgv, input%groundstate%gmaxvr )
  call generate_Gk_vectors( Gkset, kset, Gset, gkmax )
  allocate( buffer_s_alpha_matrix(n_states, first_kpt:last_kpt) )
  call obtain_s_alpha( kset, Gkset, Gset, first_kpt, buffer_s_alpha_matrix )
  call xmpi_gatherv( mpiglobal, buffer_s_alpha_matrix, s_alpha_matrix )

  print_band_character = input%groundstate%dfthalf%printBandCharacter
  if( print_band_character ) then
    l_max = min( input%groundstate%dfthalf%lmaxBandCharacter, input%groundstate%lmaxapw )
    allocate( buffer_band_character(0:l_max, natmtot, n_states_sv, first_kpt:last_kpt) )
    call obtain_band_character( kset, Gkset, first_kpt, buffer_band_character )
    buffer_pointer(1:natmtot*n_states_sv*(l_max+1), first_kpt:last_kpt) => buffer_band_character
    call xmpi_gatherv( mpiglobal, buffer_pointer, buffer_aux )
    band_character(0:l_max, 1:natmtot, 1:n_states_sv, 1:kset%nkpt) => buffer_aux
    nullify( buffer_pointer )
  end if

  if( mpiglobal%is_root ) then
    allocate( ks_eigenvalues(n_states, kset%nkpt), occupations(n_states, kset%nkpt) )
    do ik = 1, kset%nkpt
      call getevalsv( kset%vkl(:, ik), ks_eigenvalues(:, ik) )
      call getoccsv( kset%vkl(:, ik), occupations(:, ik) )
    end do
    call write_s_alpha( kset, occupations, ks_eigenvalues, s_alpha_matrix, band_character )
  end if
  if( print_band_character ) nullify( band_character )
end subroutine

subroutine set_new_occupations( is, shell, spnst, spocc, new_occupations )
  !> Species index
  integer(i32), intent(in) :: is
  !> Object that encapsulates the DFT-1/2 shell parameters
  type(dft_half_shell_parameters), intent(in) :: shell
  !> Number of states for the given species
  integer(i32), intent(in) :: spnst
  !> Original occupations for the given species
  real(dp), intent(in) :: spocc(:)
  !> Output array with the new occupations for the given species
  real(dp), intent(out) :: new_occupations(:)

  new_occupations = 0._dp
  new_occupations(1:spnst) = spocc(1:spnst)
  new_occupations(shell%number) = new_occupations(shell%number) - shell%ionization
  call terminate_if_false( all( new_occupations >= 0._dp ), &
    'Error concerning the ionization defined for species ' // to_char( is ) // &
    'A negatively occupied shell results from the present ionization condition' )
end subroutine

!> Obtain the cut-off function for DFT-1/2 potential
!> \[ f_{{\rm cut}}(r) = A \left[ 1-\left(\frac{r}{R} \right)^n \right]^{3} \]
!> where \(R\) is the cut-off radius and \(n\) is the exponent.
subroutine obtain_cut_function( radii, shell_params, cut_function )
  !> Array with the radial mesh points
  real(dp), intent(in) :: radii(:) 
  !> Object that encapsulates the DFT-1/2 shell parameters
  type(dft_half_shell_parameters), intent(in) :: shell_params
  !> Output array with the cut-off function values
  real(dp), intent(out) :: cut_function(:)

  CALL_ASSERT( size( radii ) == size( cut_function ), "radii and cut_function must have the same size" )
  where ( radii < shell_params%cut )
    cut_function = shell_params%ampl*( 1-((radii/shell_params%cut)**shell_params%exponent) )**3
  elsewhere
    cut_function = 0._dp
  end where
end subroutine

!> Obtain the \(S_\alpha\) values for all states and \(\mathbf{k}\)-points
!> \[ S_\alpha = \langle \psi_{\alpha\mathbf{k}} | V_{\rm S} | \psi_{\alpha\mathbf{k}} \rangle \]
subroutine obtain_s_alpha( kset, Gkset, Gset, first_kpt, s_alpha )
  use constants, only: zone, zzero
  use matrix_elements, only: me_init, me_mt_prepare, me_mt_mat, me_ir_mat, me_mt_alloc, me_finit
  use mod_APW_LO, only: apwfr, apword, apwordmax, lofr, lorbl, nlorb
  use mod_atoms, only: natmtot, spr, nspecies, idxas, natoms
  use mod_eigensystem, only: nmatmax
  use mod_kpointset, only: Gk_set, G_set, k_set
  use mod_muffin_tin, only: nrmt
  use mod_potential_and_density, only: rhoir, rhomt, vhalfir, vhalfmt
  use modinput, only: input
  use muffin_tin_basis, only: mt_basis_type
  !> Set of \(\mathbf{k}\)-points
  type(k_set), intent(in) :: kset
  !> Set of \(\mathbf{G}+\mathbf{k}\) vectors
  type(Gk_set), intent(in) :: Gkset
  !> Set of \(\mathbf{G}\)-points
  type(G_set), intent(in) :: Gset
  !> First \(\mathbf{k}\)-point index for this MPI process
  integer(i32), intent(in) :: first_kpt
  !> Array containing the \(S_\alpha\) values
  real(dp), intent(out) :: s_alpha(:, first_kpt:)

  integer(i32) :: i, ia, ias, ik, is, last_kpt, lmaxapw, lmaxvr, lmmaxapw, n_states, ngp
  complex(dp), allocatable :: evecfv(:, :), apwalm(:, :, :, :), vhalfig(:)
  complex(dp), allocatable :: effective_potential(:, :), mt_contribution(:, :, :)
  type(mt_basis_type) :: me_basis

  n_states = size( s_alpha, 1 )
  last_kpt = ubound( s_alpha, 2 )
  
  lmaxvr = input%groundstate%lmaxvr
  lmaxapw = input%groundstate%lmaxapw
  lmmaxapw = ( lmaxapw + 1 )**2
  allocate( apwalm(Gkset%ngkmax, apwordmax, lmmaxapw, natmtot) )

  me_basis = mt_basis_type( spr(:, 1 : nspecies), nrmt(1 : nspecies), apwfr, lofr, &
    input%groundstate%lmaxapw, apword(:, 1 : nspecies), nlorb(1 : nspecies), lorbl(:, 1 : nspecies) )
  call me_init( me_basis, lmaxvr, Gset )

  call me_mt_alloc( mt_contribution )
  do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia, is)
      ! computes gaunts times radial integrals
      call me_mt_prepare( is, ias, lmaxvr, zone, vhalfmt(:, :, ias), zzero, mt_contribution(:, :, ias) )
    end do ! natoms
  end do ! nspecies

  allocate(vhalfig(Gset%ngvec))
  call obtain_vhalfig( Gset, vhalfir, vhalfig )

  allocate( effective_potential(n_states, n_states) )
  allocate( evecfv(nmatmax, n_states) )
  do ik = first_kpt, last_kpt
    ngp = Gkset%ngk(1, ik)
    ! Get the eigenvectors from file
    call getevecfv( kset%vkl(:, ik), Gkset%vgkl(:, :, :, ik), evecfv )
    ! Matching coefficients
    call match( ngp, Gkset%gkc(:, 1, ik), Gkset%tpgkc(:, :, 1, ik), Gkset%sfacgk(:, :, 1, ik), apwalm )
    
    effective_potential = zzero
    ! mt contribution
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call me_mt_mat( is, ias, ngp, apwalm(:, :, :, ias), evecfv, zone, &
          mt_contribution(:, :, ias), zone, effective_potential )
      end do ! natoms
    end do ! nspecies
    ! ir contribution
    call me_ir_mat( Gkset, ik, evecfv, zone, vhalfig, zone, effective_potential )
    ! only diagonal terms are needed
    do i = 1, n_states
      s_alpha(i, ik) = effective_potential(i, i)%re
    end do
  end do ! ik
  call me_finit()
end subroutine

subroutine obtain_vhalfig( Gset, veffir, veffig )
  use m_zfftifc, only: zfftifc
  use mod_Gvector, only: cfunir
  use mod_kpointset, only: G_set
  !> Object that encapsulates the set of G-vectors
  type(G_set), intent(in) :: Gset
  !> Effective potential in real space (interstitial region)
  real(dp), intent(in) :: veffir(:)
  !> Effective potential in G-space
  complex(dp), intent(out) :: veffig(:)
  
  integer(i32) :: ig, ifg, ng, nr
  integer(i32), parameter :: n_cartesian = 3
  complex(dp), allocatable :: zfft(:)

  nr = size(veffir)
  ng = size(veffig)
  allocate( zfft(nr) )
  ! multiply effective potential with characteristic function
  zfft = veffir*cfunir
  ! Fourier transform to G-space
  call zfftifc( n_cartesian, Gset%ngrid, -1, zfft )
  do ig = 1, ng
    ifg = Gset%igfft(ig)
    veffig(ig) = zfft(ifg)
  end do
end subroutine

!> This subroutine computes the band character for all bands and k-points
!> It should be moved to the [[bandstructure]] module in the future
!> But there are no tests for that module yet
subroutine obtain_band_character( kset, Gkset, first_kpt, bc )
  use mod_APW_LO, only: apwordmax
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_kpointset, only: Gk_set, k_set
  use mod_muffin_tin, only: lmmaxapw
  !> Set of \(\mathbf{k}\)-points
  type(k_set), intent(in) :: kset
  !> Set of \(\mathbf{G}+\mathbf{k}\) vectors 
  type(Gk_set), intent(in) :: Gkset
  !> First \(\mathbf{k}\)-point index for this MPI process
  integer(i32), intent(in) :: first_kpt
  !> Array containing the band characters
  real(dp), contiguous, intent(out) :: bc(0:, :, :, first_kpt:)

  integer(i32) :: ik, last_kpt, n_atoms, n_states, n_states_sv, ngp
  complex(dp), allocatable :: evecfv(:, :), evecsv(:, :), apwalm(:, :, :, :)

  n_states = nstfv
  n_states_sv = size( bc, 3 )
  last_kpt = ubound( bc, 4 )
  allocate( evecfv(nmatmax, n_states), evecsv(n_states, n_states_sv) ) 
  n_atoms = size( bc, 2 )
  allocate( apwalm(Gkset%ngkmax, apwordmax, lmmaxapw, n_atoms) )
  do ik = first_kpt, last_kpt
    ! Get the eigenvectors from file
    call getevecfv( kset%vkl(:, ik), Gkset%vgkl(:, :, :, ik), evecfv )
    call getevecsv( kset%vkl(:, ik), Gkset%vgkl(:, :, :, ik), evecsv )
    ! Matching coefficients
    ngp = Gkset%ngk(1, ik)
    call match( ngp, Gkset%gkc(:, 1, ik), Gkset%tpgkc(:, :, 1, ik), Gkset%sfacgk(:, :, 1, ik), apwalm )
    call get_band_character_ik( ngp, evecfv, evecsv, apwalm, bc(:, :, :, ik) )
  end do
end subroutine

!> This subroutine computes the band character for all bands and a given k-point
subroutine get_band_character_ik( ngk, evecfv, evecsv, apwalm, bc_ik )
  use mod_atoms, only: nspecies, idxas, natoms
  use mod_spin, only: nspinor
  !> Number of \(\mathbf{G}+\mathbf{k}\) vectors
  integer(i32), intent(in) :: ngk
  !> KS eigenvectors for the given \(\mathbf{k}\)-point (in LAPW+LO basis)
  complex(dp), contiguous, intent(in) :: evecfv(:, :)
  !> KS eigenvectors for the given \(\mathbf{k}\)-point (in second-variational basis)
  complex(dp), contiguous, intent(in) :: evecsv(:, :)
  !> Matching coefficients for the given \(\mathbf{k}\)-point
  complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :)
  !> Output array with the band character for the given \(\mathbf{k}\)-point
  real(dp), intent(out) :: bc_ik(0:, :, :)
  
  integer(i32) :: ia, ias, is, ispn, ist, l, lm, l_max, lm_max, m, n_states
  real(dp) :: aux
  complex(dp), allocatable :: dmat(:, :, :, :, :)

  l_max = ubound( bc_ik, 1 )
  lm_max = (l_max+1)**2
  n_states = size( bc_ik, 3 )
  allocate( dmat(lm_max, lm_max, nspinor, nspinor, n_states) )
  ! average band character over spin and m for all atoms
  do is = 1, nspecies
    do ia = 1, natoms (is)
      ias = idxas(ia, is)
      ! generate the diagonal of the density matrix
      call gendmat(.true., .true., 0, l_max, is, ia, ngk, apwalm, evecfv, evecsv, lm_max, dmat)
      do ist = 1, n_states
        lm = 1
        do l = 0, l_max
          aux = 0._dp
          do m = -l, l
            do ispn = 1, nspinor
              aux = aux + dmat(lm, lm, ispn, ispn, ist)%re
            end do
            lm = lm + 1
          end do
          bc_ik(l, ias, ist) = aux
        end do
      end do
    end do
  end do
end subroutine

subroutine write_VS_to_file( is, shell_params_array, radii, Vatom, Vion, VS )
  !> Species index
  integer(i32), intent(in) :: is
  !> Object array that encapsulates the DFT-1/2 shell parameters
  type(dft_half_shell_parameters), intent(in) :: shell_params_array(:)
  real(dp), intent(in) :: radii(:), Vatom(:), Vion(:), VS(:)

  integer(i32) :: i_unit
  character(len=20) :: file_name

  write( file_name, '("VS_S", I2.2, ".OUT")') is
  open( newunit=i_unit, file=trim(file_name), action='WRITE', form='FORMATTED')
  call write_VS_file_header( is, i_unit, shell_params_array )
  call write_radii_and_potentials( i_unit, radii, Vatom, Vion, VS )
  close( i_unit )
end subroutine

subroutine write_VS_file_header( is, i_unit, shell_params_array )
  !> Species index
  integer(i32), intent(in) :: is
  !> File unit to write the output
  integer(i32), intent(in) :: i_unit
  !> Object array that encapsulates the DFT-1/2 shell parameters
  type(dft_half_shell_parameters), intent(in) :: shell_params_array(:)

  integer(i32) :: i, n_shell
  n_shell = size( shell_params_array )

  write( i_unit, "(A, I4)") "Species: ", is
  write( i_unit, "(I4, A)") n_shell, " shell(s) must be ionized. Its/Their parameters(s) are listed below."
  write( i_unit, "(A6, 4A15)") "shell", "ionization", "cut_amplitude", "cut_radius", "cut_exponent"
  associate( shell => shell_params_array )
    do i = 1, n_shell
      write( i_unit,'(I6, 3F15.6, I15)') shell(i)%number, shell(i)%ionization, &
        shell(i)%ampl, shell(i)%cut, shell(i)%exponent
    end do
  end associate
end subroutine

subroutine write_radii_and_potentials( i_unit, radii, Vatom, Vion, VS )
  integer(i32), intent(in) :: i_unit
  real(dp), intent(in) :: radii(:), Vatom(:), Vion(:), VS(:)

  integer(i32) :: ir

  write( i_unit, "(A, I8)") "number of radial points: ", size( radii )
  write( i_unit, "(4A20)") "r", "Vatom", "Vion", "VS"
  do ir = 1, size( radii )
    write( i_unit,"(4ES20.6E2)") radii(ir), Vatom(ir), Vion(ir), VS(ir)
  end do
end subroutine

subroutine write_s_alpha( kset, occupations, ks_eigenvalues, s_alpha, band_character )
  use mod_kpointset, only: k_set
  use to_char_conversion, only: to_char
  type(k_set), intent(in) :: kset
  real(dp), intent(in) :: occupations(:, :)
  real(dp), intent(in) :: ks_eigenvalues(:, :)
  real(dp), intent(in) :: s_alpha(:, :)
  real(dp), intent(in), optional, contiguous :: band_character(0:, :, :, :)

  character(len=:), allocatable :: buffer, character_string, header
  character(len=*), parameter :: output_name = "DFT_HALF_NSCF.OUT"
  character(len=*), parameter :: string_format_header = '(I4, 3G13.5, A)'
  character(len=*), parameter :: string_format = '(I4, 3G13.5)'
  integer(i32), parameter :: string_len = 4 + 3 * 13  ! I4 + 3 * G13.5
  integer(i32) :: char_len, funit, i, ik, j, line_len, n_states, n_atoms, n_characters, pos
  logical :: print_band_character

  n_states = size( s_alpha, 1 )
  print_band_character = present( band_character )
  header = "(state, eigenvalue,   occupancy,   Salpha)        "

  line_len = string_len
  if( print_band_character ) then
    n_atoms = size( band_character, 2 )
    n_characters = ubound( band_character, 1 )
    do i = 1, n_atoms
      header = header // ", at" // to_char(i) //":"
      do j = 0, n_characters
        header = header // "l=" // to_char(j) // " "
      end do
    end do
    character_string = "(" // to_char(n_characters+1) // "F8.4)"
    char_len = 8 * (n_characters + 1)  ! "F8.4" per character, n_characters+1 of them
    line_len = line_len + n_atoms * char_len
  end if
  allocate( character(len=line_len) :: buffer )
  open( newunit=funit, file=output_name, status='replace', action='write' )
  write( funit, '(I6, A)') kset%nkpt, " : nkpt"
  write( funit, '(I6, A)') n_states, " : nstsv"
  do ik = 1, kset%nkpt
    write( funit, string_format_header ) ik, kset%vkl(:, ik), " : k-point, vkl"
    write( funit, '(A)') header
    do i = 1, n_states
      pos = string_len
      write( buffer(1:pos), string_format ) i, ks_eigenvalues(i, ik), occupations(i, ik), s_alpha(i, ik)
      if( print_band_character ) then
        do j = 1, n_atoms
          write( buffer(pos+1:pos+char_len), character_string ) band_character(:, j, i, ik)
          pos = pos + char_len
        end do
      end if
      write( funit, '(a)' ) trim(buffer)
    end do
    write( funit, * ) ""
  end do
  close( funit )
end subroutine
end module
