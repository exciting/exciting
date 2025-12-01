! TODO(Ronaldo): Refactor to reduce the number of global variables
!> Module that manages the hamiltonian matrix in RT-TDDFT
module rttddft_Hamiltonian
  use asserts, only: assert
  use constants, only: fourpi, y00, zi, zone, zzero
  use matrix_elements, only: me_mt_alloc, me_mt_prepare, me_mt_mat, me_ir_alloc, me_ir_prepare, me_ir_mat
  use mod_atoms, only: atposc, idxas, natoms, natmtot, nspecies
  use mod_gvector, only: cfunig
  use mod_kpointset, only: Gk_set
  use mod_lattice, only: omega
  use mod_muffin_tin, only: rmt
  use mod_potential_and_density, only: meffig, veffig, veffir, veffmt
  use modinput, only: input
  use physical_constants, only: alpha, c
  use precision, only: dp, i32
  use rttddft_arrays, only: generic_matrix_set, hermitian_matrix_set
  use rttddft_Overlap, only: overlap_set
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_hamiltonian, &
    Timing_Ehrenfest, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use xlapack, only: scaled_add

  implicit none

  private
  
  integer(i32), parameter :: n_cartesian = 3
  real(dp), parameter :: eps_scissor = 1.e-10_dp

  !> Type to encapsulate the set of hamiltonian matrices (which can be not hermitian)
  type, public :: hamiltonian_set
    private
    !> If `.true.`, use the IPA
    logical, private :: IPA = .false.
    !> \(\mathbf{k}\)-dependent dimensions array
    integer(i32), allocatable, public :: dims(:)
    !> Hamiltonian matrix at time \(t = 0\)
    type(hermitian_matrix_set), public :: H_0
    !> Hamiltonian matrix at time \(t\)
    class(generic_matrix_set), allocatable, public :: H_t
    !> Hamiltonian matrix at previous time \(t - \Delta t \)
    class(generic_matrix_set), allocatable, public :: H_t_minus_dt
    !> Matrix elements of the KS potential at time \(t = 0\)
    type(hermitian_matrix_set), public :: V_KS_0
    !> Hamiltonian eigenvalues at time \(t = 0\)
    real(dp), public, allocatable :: initial_eigenvalues(:, :)
    !> If `.true.`, explicit evaluation should be performed in [[hamiltonian_set_calculate]]
    logical, private :: explicit_evaluation_needed = .true.
    !> `mathcalH` gives the impact of an ion displacement on the hamiltonian matrix
    !> \[ \left[ \left\langle 
    !> \frac{\partial \phi_{\mu'}^{\mathbf{k}}}{\partial \mathbf{R}_J}
    !> \Bigg|\hat{H}\Bigg|\phi_{\mu}^{\mathbf{k}}\right\rangle +
    !> \left\langle\phi_{\mu'}^{\mathbf{k}}\Bigg|\hat{H}\Bigg|\frac{\partial 
    !> \phi_{\mu}^{\mathbf{k}}}{\partial \mathbf{R}_J}\right\rangle \right] 
    !> \]
    complex(dp), public, allocatable  :: mathcalH(:, :, :, :, :)
  contains
    private
    procedure, private :: add_external_coupling_berry_phase => add_external_coupling_berry_phase
    procedure, private :: add_external_coupling_velocity_gauge => add_external_coupling_velocity_gauge
    generic, public    :: add_external_coupling => add_external_coupling_berry_phase, &
                                                   add_external_coupling_velocity_gauge
    procedure, public  :: allocate => hamiltonian_set_allocate
    procedure, private :: calculate_in_lapw_basis => hamiltonian_set_calculate_lapw_basis
    procedure, private :: calculate_in_ks_basis => hamiltonian_set_calculate_ks_basis
    procedure, private :: calculate_ik => hamiltonian_set_calculate_lapw_basis_ik
    procedure, public  :: calculate => hamiltonian_set_calculate
    procedure, public  :: copy_H_t => hamiltonian_set_copy_H_t
    procedure, public  :: adjust_eigenvalues_with_scissor_shift => hamiltonian_set_adjust_eigenvalues_with_scissor_shift
    final              :: destructor
  end type

contains

  subroutine hamiltonian_set_allocate( this, max_dimension, ki, dims, n_states, &
      allocate_H_past, evolve_H0, is_LAPWLO_basis, MD, is_IPA )
    class(hamiltonian_set), intent(inout) :: this
    !> Maximum dimension of the matrices for all \( \mathbf{k} \)-points
    integer(i32), intent(in) :: max_dimension
    !> First \( \mathbf{k} \)-point managed by MPI process
    integer(i32), intent(in) :: ki
    !> \( \mathbf{k} \)-dependent dimensions
    integer(i32), intent(in) :: dims(ki:)
    !> Number of KS states (n_occupied + n_empty). Relevant only with `is_LAPWLO_basis` = `.False.`
    integer(i32), intent(in) :: n_states
    !> if `.True`, allocate `H_t_minus_dt`
    logical, intent(in) :: allocate_H_past
    !> if `.True`, evolve explicitly time-independent `H_0`
    logical, intent(in) :: evolve_H0
    !> if `.True`, LAPW+lo basis is used
    logical, intent(in) :: is_LAPWLO_basis
    !> if `.True`, MD is carried out
    logical, intent(in) :: MD ! TODO: two logicals should actually be passed: molecular_dynamics%on and molecular_dynamics%valence_corrections
    !> if `.True`, independent particle approximation is employed
    logical, intent(in) :: is_IPA

    logical :: is_KS_basis
    integer(i32) :: kf

    call assert( max_dimension >= maxval( dims ), 'm must be >= maxval( dims )' )
    this%IPA = is_IPA
    is_KS_basis = .not. is_LAPWLO_basis
    if ( is_KS_basis ) call assert ( max_dimension == n_states, 'max_dimension must be equal to n_states' )
    kf = ubound( dims, 1 )
    if( MD ) then ! H_t and H_t_minus_dt matrices are not hermitian
      allocate( generic_matrix_set :: this%H_t )
      allocate( generic_matrix_set :: this%H_t_minus_dt )
    else
      allocate( hermitian_matrix_set :: this%H_t )
      allocate( hermitian_matrix_set :: this%H_t_minus_dt )
    end if
    call this%H_t%allocate_array( [1, 1, ki], [max_dimension, max_dimension, kf] )
    if( allocate_H_past ) call this%H_t_minus_dt%allocate_array( [1, 1, ki], [max_dimension, max_dimension, kf] )
    if( ( .not. evolve_H0 ) .and. ( is_LAPWLO_basis ) ) call this%H_0%allocate_array( [1, 1, ki], [max_dimension, max_dimension, kf] ) 
    if( MD ) allocate( this%mathcalH(max_dimension, max_dimension, n_cartesian, natmtot, ki:kf) )
    if ( is_KS_basis ) call this%V_KS_0%allocate_array( [1, 1, ki], [max_dimension, max_dimension, kf] )
    if( allocated( this%dims ) ) deallocate( this%dims )
    if( is_KS_basis ) then 
      allocate( this%dims(ki:kf), source = max_dimension ) 
    else ! This represents the LAPW basis
      allocate( this%dims(ki:kf), source = dims(ki:kf) )
    end if
    if( allocated( this%initial_eigenvalues ) ) deallocate( this%initial_eigenvalues )
    allocate( this%initial_eigenvalues(n_states, ki:kf) )
  end subroutine

  !> Copy `H_t` into `H_t_minus_dt`
  subroutine hamiltonian_set_copy_H_t( this )
    class(hamiltonian_set), intent(inout) :: this
    call this%H_t_minus_dt%copy_from( this%H_t )
  end subroutine

  subroutine destructor( this )
    type(hamiltonian_set), intent(inout) :: this

    if( allocated( this%H_t ) ) call this%H_t%deallocate_if_allocated()
    if( allocated( this%H_t_minus_dt ) ) call this%H_t_minus_dt%deallocate_if_allocated()
    call this%H_0%deallocate_if_allocated()
    call this%V_KS_0%deallocate_if_allocated()
    if( allocated( this%initial_eigenvalues ) ) deallocate( this%initial_eigenvalues )
    if( allocated( this%mathcalH ) ) deallocate( this%mathcalH )
    if( allocated( this%dims ) ) deallocate( this%dims )
  end subroutine

  !> Interface to decide if calculate in LAPW or KS basis and call the corresponding subroutines.
  subroutine hamiltonian_set_calculate( this, l_max_pot, apwalm, Gkset, ks_lapwlo_transition_matrix,&
    printTimings, t_ham, a_tot, obtain_mathcalH, t_MD )
    class(hamiltonian_set), intent(inout) :: this
    !> Maximal value of l in spherical harmonics expansion of DFT potential
    integer(i32), intent(in) :: l_max_pot
    !> Matching coefficients of the (L)APWs (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> KS-LAPW+lo transition matrix (nmatmax, n_basis_ks, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: ks_lapwlo_transition_matrix(:, :, :)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(inout) :: t_ham
    !> Total vector potential
    type(Vector_Potential_Field), optional, intent(in) :: a_tot
    !> If `.true.`, obtain `mathcalH`
    logical, optional, intent(in) :: obtain_mathcalH
    !> Object that packs information about timings related to MD
    type(Timing_Ehrenfest), optional, intent(out) :: t_MD

    if( present( ks_lapwlo_transition_matrix ) ) then
      call this%calculate_in_ks_basis( l_max_pot, apwalm, Gkset, &
        ks_lapwlo_transition_matrix, printTimings, t_ham )
    else  
      call this%calculate_in_lapw_basis( l_max_pot, apwalm, Gkset, &
        printTimings, t_ham, a_tot, obtain_mathcalH, t_MD )
    end if
    ! with IPA, the Hamiltonian should only be calculated from scratch at t = 0
    if ( this%IPA ) this%explicit_evaluation_needed = .false.
  end subroutine

  !> Obtain the field-free hamiltonian at time \( t \) in the LAPW+lo basis.
  subroutine hamiltonian_set_calculate_lapw_basis( this, l_max_pot, apwalm, Gkset, &
    printTimings, t_ham, a_tot, obtain_mathcalH, t_MD )
    class(hamiltonian_set), intent(inout) :: this
    !> Maximal value of l in spherical harmonics expansion of DFT potential
    integer(i32), intent(in) :: l_max_pot
    !> Matching coefficients of the (L)APWs (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(inout) :: t_ham
    !> Total vector potential
    type(Vector_Potential_Field), optional, intent(in) :: a_tot
    !> If `.true.`, obtain `mathcalH`
    logical, optional, intent(in) :: obtain_mathcalH
    !> Object that packs information about timings related to MD
    type(Timing_Ehrenfest), optional, intent(out) :: t_MD

    integer(i32) :: ik, first_kpt, last_kpt, ias, ia, is, n_MT_radial_points
    real(dp) :: t_i, t_f, t_aux
    logical :: timings_general, timings_detailed, valence_relativity
    integer(i32), parameter :: l_max_kin = 0 ! Overlap operator is equal to (1.0/y00)*Y_{00}, it has only l=0 component
    integer(i32), parameter :: lm_max = (l_max_kin+1)**2
    real(dp), parameter :: aux_factor = 0.5_dp*alpha**2*y00
    real(dp), parameter :: kin_operator_00 = (0.5_dp / y00)
    real(dp), allocatable :: kin_mt(:, :)
    complex(dp), allocatable :: mt_part(:, :, :)

    first_kpt = lbound( this%H_t%array, 3)
    last_kpt = ubound( this%H_t%array, 3 )
    ! Check optional arguments
    timings_general = .False.
    timings_detailed = .False.
    if( present( printTimings ) ) call printTimings%get( timings_general, timings_detailed )

    ! sanity checks
    if( timings_general ) call assert( present( t_ham ) .or. present( t_MD ), &
      't_ham or t_MD must be present when general timing is desired' )
    if( timings_detailed ) call assert( timings_general, 'timings_general must be true if timings_detailed is true')

    if( timings_general ) call timesec( t_i )
    if( this%explicit_evaluation_needed ) then
      ! Interface to input variables
      valence_relativity = ( trim( input%groundstate%ValenceRelativity ) /= 'none' )

      ! MT part (prepare)
      call me_mt_alloc( mt_part )
      n_MT_radial_points = size( veffmt, 2 )
      allocate( kin_mt(lm_max, n_MT_radial_points), source = kin_operator_00 ) ! Kinetic operator = (0.5/y00)*Y_{00}
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! If valence_relativity is true, then kin_mt = (0.5/y00) / (1.0 - 0.5 * alpha^2 * veffmt * y00)
          if( valence_relativity ) &
            kin_mt = kin_operator_00 / (1.0_dp - aux_factor * veffmt(1:lm_max, :, ias) )
          call me_mt_prepare( is, ias, l_max_kin, zone, kin_mt, zone, mt_part(:, :, ias), gradient_product=.true. )
          call me_mt_prepare( is, ias, l_max_pot, zone, veffmt(:, :, ias), zone, mt_part(:, :, ias) )
        end do
      end do
      if( timings_detailed .and. present( t_ham ) ) then
        t_aux = t_i
        call timesec_RTTDDFT( t_aux, t_ham%hmlint )
      end if

      if( valence_relativity ) then 
        do ik = first_kpt, last_kpt
          call this%calculate_ik( ik, Gkset, apwalm(:, :, :, :, ik-first_kpt+1), &
            mt_part, meffig, a_tot, obtain_mathcalH )
        end do
      else
        do ik = first_kpt, last_kpt
          call this%calculate_ik( ik, Gkset, apwalm(:, :, :, :, ik-first_kpt+1), &
            mt_part, cfunig, a_tot, obtain_mathcalH )
        end do
      end if
    else
      call this%H_0%assert_allocated()
      call this%H_t%copy_from( this%H_0 )
    end if

    if( timings_general ) then
      call timesec( t_f )
      if( present( t_ham ) ) t_ham%total = t_f - t_i
      if( timings_detailed .and. present( t_MD ) ) t_MD%ham = t_f - t_i
    end if
  end subroutine hamiltonian_set_calculate_lapw_basis

  !> Calculate the hamiltonian matrix in the LAPW+lo basis for a given \( \mathbf{k} \)-point.
  subroutine hamiltonian_set_calculate_lapw_basis_ik( this, ik, Gkset, apwalm_ik, mt_part, &
      kin_ir, a_tot, obtain_mathcalH )
    class(hamiltonian_set), intent(inout) :: this
    !> ik: the index of the k-point considered
    integer(i32), intent(in) :: ik
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Matching coefficients of the (L)APWs at the current k-point (ngkmax, apwordmax, lmmaxapw, natmtot)
    complex(dp), contiguous, intent(in) :: apwalm_ik(:, :, :, :)
    !> MT part of the overlap matrix
    complex(dp), contiguous, intent(in) :: mt_part(:, :, :)
    !> IR - kinetic part
    complex(dp), contiguous, intent(in) :: kin_ir(:)
    !> `x`, `y`, and `z` components of the (total) vector potential
    !> Total vector potential
    type(Vector_Potential_Field), optional, intent(in) :: a_tot
    !> If `.true.`, obtain `mathcalH`
    logical, optional, intent(in) :: obtain_mathcalH

    integer(i32) :: ias, ia, is
    complex(dp), allocatable :: tmp(:, :)
    logical :: get_mathcalH

    get_mathcalH = .false.
    if( present( obtain_mathcalH ) ) get_mathcalH = obtain_mathcalH .and. allocated( this%mathcalH )
    this%H_t%array(:, :, ik) = zzero
    ! It is better to split the two cases (even with some code duplication): 
    ! - to avoid an "if" inside the double loop over atoms
    ! - to avoid allocating "tmp" when not needed (this can be a large array)
    associate( np => Gkset%ngk(1, ik) )
    if( get_mathcalH ) then
      call assert( present( a_tot ), "a_tot must be present" )
      tmp = this%H_t%array(:, :, ik)
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, np, apwalm_ik(:, :, :, ias), zone, &
            mt_part(:, :, ias), zzero, tmp )
          this%H_t%array(:, :, ik) = this%H_t%array(:, :, ik) + tmp
          call update_mathcalH_ik_ias( is, ia, np, Gkset%vgkc(:, :, 1, ik), tmp, &
            a_tot%components, this%mathcalH(:, :, :, ias, ik) )
        end do
      end do
    else
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, np, apwalm_ik(:, :, :, ias), zone, &
            mt_part(:, :, ias), zone, this%H_t%array(:, :, ik) )
        end do
      end do
    end if
    end associate
    call me_ir_mat( Gkset, ik, zone, veffig, zone, this%H_t%array(:, :, ik) )
    call me_ir_mat( Gkset, ik, zone/2, kin_ir, zone, this%H_t%array(:, :, ik), gradient_product=.true. )
  end subroutine

  !> (private) Update the `mathcalH` matrix for a given atom and \( \mathbf{k} \)-point.
  subroutine update_mathcalH_ik_ias( is, ia, n_pw, gplusk_cart, ham_ik_ias, a_tot, mathcalH_ik_ias )
    !> Species index
    integer(i32), intent(in) :: is
    !> Atom index
    integer(i32), intent(in) :: ia
    !> Number of plane waves
    integer(i32), intent(in) :: n_pw
    !> G+k vector in Cartesian coordinates (3, ngp)
    real(dp), contiguous, intent(in) :: gplusk_cart(:, :)
    !> MT part of hamiltonian for a given atom and k-point
    complex(dp), contiguous, intent(in) :: ham_ik_ias(:, :)
    !> `x`, `y`, and `z` components of the (total) vector potential
    real(dp), optional, intent(in) :: a_tot(n_cartesian)
    !> `mathcalH` matrix, given an atom and k-point 
    complex(dp), intent(inout) :: mathcalH_ik_ias(:, :, :)

    integer(i32) :: ig, igl, i_cart
    real(dp) :: aux, diff(n_cartesian), a_scaled(n_cartesian), t1, t2, t3, t4
    complex(dp) :: t5

    associate( m => size(mathcalH_ik_ias, 1) )
      call assert( size(gplusk_cart, 1) == n_cartesian, "gplusk_cart must have n_cartesian elements along 1st dim" )
      call assert( size(gplusk_cart, 2) >= n_pw, "gplusk_cart must at least n_pw elements along 2nd dim" )
      call assert( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "mathcalH_ik_ias must have shape [m, m, n_cartesian]" )
      call assert( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "p_MT_ias_ik must have shape [m, m, n_cartesian]" )
      call assert( n_pw <= m, "n_pw must be <= m" )
      a_scaled = a_tot / c
      t1 = fourpi*(rmt(is)**3)/omega
      do i_cart = 1, n_cartesian
        do ig = 1, n_pw
          do igl = 1, n_pw
            diff = gplusk_cart(:, ig) - gplusk_cart(:, igl)
            aux = diff(i_cart)
            t2 = dot_product(0.5_dp*gplusk_cart(:, igl) + a_scaled, gplusk_cart(:, ig) )
            t3 = rmt(is)*sqrt( diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3) )
            if( igl /= ig ) t3 = ( sin(t3)-t3*cos(t3) )/( t3**3 ) ! for igl == ig, t3 = 0
            t4 = dot_product( diff, atposc(:, ia, is) )
            t5 = cmplx( cos(t4), sin(t4), kind(dp) )
            mathcalH_ik_ias(igl, ig, i_cart) = mathcalH_ik_ias(igl, ig, i_cart) + &
              zi*aux*( ham_ik_ias(igl, ig) - t1*t2*t3*t5 )
          end do
          aux = gplusk_cart(i_cart, ig)
          do igl = n_pw + 1, m
            mathcalH_ik_ias(igl, ig, i_cart) = mathcalH_ik_ias(igl, ig, i_cart) + &
              zi*aux*ham_ik_ias(igl, ig)
          end do
        end do
        do ig = n_pw + 1, m
          do igl = 1, n_pw
            aux = -gplusk_cart(i_cart, igl)
            mathcalH_ik_ias(igl, ig, i_cart) = mathcalH_ik_ias(igl, ig, i_cart) + &
              zi*aux*ham_ik_ias(igl, ig) 
          end do
        end do
      end do
    end associate
  end subroutine

  !> Obtain the explicitly field-independent 
  !> Hamiltonian \( H(t) \) in the KS basis as follows:
  !> \[
  !> H(t) = H_{\rm init} + V_{\rm eff}(t) - V_{\rm eff}(t = 0), 
  !> \]
  !> where
  !> \[ 
  !> H_{\rm init} = \frac{{\bf p}^2}{2} + V_{\rm nucl} + V_{\rm eff}(t = 0)
  !> \]
  !> is time-independent. \( V_{\rm eff}(t) \) is time-dependent through the 
  !> time-dependent charge density.
  subroutine hamiltonian_set_calculate_ks_basis( this, lmaxvr, apwalm, Gkset, &
      ks_lapwlo_transition_matrix, printTimings, t_ham )
    class(hamiltonian_set), intent(inout) :: this
    !> Maximal value of l in spherical harmonics expansion of DFT potential
    integer(i32), intent(in) :: lmaxvr
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> Set of \( \mathbf{G} + \mathbf{k} \) vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> KS-LAPW+lo transition matrix (nmatmax, n_basis_ks, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: ks_lapwlo_transition_matrix(:, :, :)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(out) :: t_ham

    integer(i32) :: ik, first_kpt, last_kpt, shift, n_basis, is, ia, ias, ngp
    complex (dp), allocatable :: local_effective_potential(:, :), mt_contribution(:, :, :)
    logical :: timings_general, timings_detailed
    real(dp) :: ti

    first_kpt = lbound( this%H_t%array, 3 )
    last_kpt = ubound( this%H_t%array, 3 )
    n_basis = size( this%H_t%array, 1 )
    shift = first_kpt - 1

    timings_general = .False.
    timings_detailed = .False.
    if ( present( printTimings ) ) call printTimings%get( timings_general, timings_detailed )
    if( timings_general ) call assert( present( t_ham ), &
      't_ham must be present when general timing is desired' )
    if( timings_detailed ) call assert( timings_general, 'timings_general must be true if timings_detailed is true')

    if( timings_general ) call timesec( ti )

    this%H_t%array = zzero
    call this%H_t%copy_from( this%initial_eigenvalues )

    if( this%explicit_evaluation_needed ) then
      call this%H_t%subtract( this%V_KS_0 )

      call me_mt_alloc( mt_contribution )
      allocate( local_effective_potential(n_basis, n_basis) )
    
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! computes gaunts times radial integrals
          call me_mt_prepare( is, ias, lmaxvr, zone, veffmt(:, :, ias), zzero, &
            mt_contribution(:, :, ias) )
        end do ! natoms
      end do ! nspecies

      do ik = first_kpt, last_kpt
        ngp = Gkset%ngk(1, ik)
        local_effective_potential = zzero
        ! mt contribution
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            call me_mt_mat( is, ias, ngp, apwalm(:, :, :, ias, ik-shift), &
              ks_lapwlo_transition_matrix(:, :, ik-shift), zone, mt_contribution(:, :, ias), &
              zone, local_effective_potential )
          end do ! natoms
        end do ! nspecies
        ! ir contribution
        call me_ir_mat( Gkset, ik, ks_lapwlo_transition_matrix(:, :, ik-shift), &
          zone, veffig, zone, local_effective_potential )
        ! final hamiltonian
        this%H_t%array(:, :, ik) = this%H_t%array(:, :, ik) + local_effective_potential
      end do ! ik
    end if ! this%explicit_evaluation_needed

    if( timings_general ) call timesec_RTTDDFT( ti, t_ham%total )
  end subroutine
  
  !> Add the pre-calculated length gauge interaction term to the Hamiltonian
  subroutine add_external_coupling_berry_phase( this, external_coupling_length_gauge )
    class(hamiltonian_set), intent(inout) :: this
    !> Length gauge interaction matrix (n_basis, n_basis, n_kpts)
    complex(dp), contiguous, intent(in) :: external_coupling_length_gauge(:, :, :)

    call assert( all( shape( external_coupling_length_gauge ) == shape( this%H_t%array ) ), &
      'external_coupling_length_gauge and hamiltonian have incompatible dimensions' )

    this%H_t%array = this%H_t%array + external_coupling_length_gauge
  end subroutine

  !> Add the velocity gauge interaction term \( {\bf p} \cdot {\bf A}(t) / c \) to 
  !> the Hamiltonian at time \( t \).
  ! TODO: is the space-uniform A^2 term needed here?
  subroutine add_external_coupling_velocity_gauge( this, a_tot, overlap, pmat )
    class(hamiltonian_set), intent(inout) :: this
    !> Total vector potential
    type(Vector_Potential_Field), intent(in) :: a_tot
    !> Overlap matrix
    class(overlap_set), intent(in) :: overlap
    !> Momentum matrix elements (n_basis, n_basis, 3, n_kpts)
    complex(dp), contiguous, intent(in) :: pmat(:, :, :, :)

    real(dp), parameter :: interaction_tol = 1.e-14_dp
    integer(i32) :: ik, i
    real(dp) :: a_scaled(n_cartesian), fact

    a_scaled = a_tot%components / c
    fact = 0.5_dp * dot_product( a_scaled, a_scaled )
    if ( fact < interaction_tol ) return
    associate( m => size( this%H_t%array, 1 ), n_kpts => size( this%H_t%array, 3 ) )
      call assert( size( pmat, 4 ) == n_kpts, "pmat and hamiltonian have different n_kpts" )
      if( overlap%is_identity() ) then
        ! TODO: use DO CONCURRENT here after ifort2021 support is dropped
        do ik = lbound(this%H_t%array, 3), ubound(this%H_t%array, 3)
          do i = 1, m
            this%H_t%array(i, i, ik) = this%H_t%array(i, i, ik) + cmplx( fact, kind = dp )
          end do
        end do
      else
        call assert( all( shape( overlap%array ) == shape( this%H_t%array ) ), &
          "overlap and hamiltonian must have same shape" )
        call scaled_add( fact, overlap%array, this%H_t%array )
      end if
      do i = 1, n_cartesian
        call scaled_add( a_scaled(i), pmat(:, :, i, :), this%H_t%array )
      end do
    end associate
  end subroutine

  !> Adjust the initial_eigenvalues arrays, shifting the conduction band energies upwards by \( \Delta E \).
  subroutine hamiltonian_set_adjust_eigenvalues_with_scissor_shift( this, scissor_shift, first_unoccupied )
    class(hamiltonian_set), intent(inout) :: this
    !> Energy shift  \( \Delta E \) for scissor operator
    real(dp), intent(in) :: scissor_shift
    !> Position of the first unoccupied state
    integer, intent(in) :: first_unoccupied

    if ( scissor_shift > eps_scissor ) &
      this%initial_eigenvalues(first_unoccupied:, :) = this%initial_eigenvalues(first_unoccupied:, :) + scissor_shift
  end subroutine

end module rttddft_Hamiltonian
