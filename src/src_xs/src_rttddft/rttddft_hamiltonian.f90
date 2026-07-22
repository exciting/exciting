! TODO(Ronaldo): Refactor to reduce the number of global variables
!> Module that manages the hamiltonian matrix in RT-TDDFT
module rttddft_Hamiltonian
#include "asserts.fpp"
  use constants, only: fourpi, real_zero, y00, zi, zone, zzero
  use generation_wavefunction, only: generate_basisfunction_secondvariation_MT
  use math_utils, only: is_hermitian
  use matrix_elements, only: me_mt_alloc, me_mt_prepare, me_mt_mat, me_ir_alloc, me_ir_prepare, me_ir_mat
  use mod_atoms, only: atposc, idxas, natoms, natmtot, nspecies, spr
  use mod_gvector, only: cfunig
  use mod_kpointset, only: Gk_set
  use mod_lattice, only: omega
  use mod_muffin_tin, only: nrcmtmax, nrmt, rcmt, rmt
  use mod_potential_and_density, only: meffig, veffig, veffir, veffmt
  use mod_variation, only: variation_multiplication
  use modinput, only: input, isspinorb
  use physical_constants, only: alpha, c
  use precision, only: dp, i32
  use rttddft_arrays, only: generic_matrix_set, hermitian_matrix_set
  use rttddft_Overlap, only: overlap_set
  use rttddft_pmat, only: pmat_set
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_hamiltonian, &
    Timing_Ehrenfest, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use to_char_conversion, only: to_char
  use xlapack, only: matrix_multiply, scaled_add

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
    !> Matrix of the scissor shift operator
    type(hermitian_matrix_set), private :: scissor_matrix
    !> Hamiltonian eigenvalues at time \(t = 0\)
    real(dp), public, allocatable :: initial_eigenvalues(:, :)
    !> If `.true.`, explicit evaluation should be performed in [[hamiltonian_set_calculate]]
    logical, private :: explicit_evaluation_needed = .true.
    !> If `.true.`, the SOC part of the hamiltonian should be recalculated every time H is calculated
    logical, private :: obtainSOC = .true.
    !> If `.true.`, Hamiltonian is built in the LAPW+lo basis set
    logical, private :: lapwlo_basis = .true.
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
    procedure, public  :: build_lapwlo_scissor_matrix => hamiltonian_set_build_lapwlo_scissor_matrix
    procedure, public  :: adjust_eigenvalues_with_scissor_shift => hamiltonian_set_adjust_eigenvalues_with_scissor_shift
    procedure, public  :: represented_in_lapwlo => hamiltonian_set_represented_in_lapwlo
    final              :: destructor
  end type

contains

  subroutine hamiltonian_set_allocate( this, max_dimension, ki, dims, n_states, &
      allocate_H_past, evolve_H0, is_LAPWLO_basis, MD, is_IPA, scissor_shift, update_SOC )
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
    !> Requested scissor energy correction value
    real(dp), intent(in) :: scissor_shift
    !> if `.True`, SOC should be updated every time H is calculated
    logical, intent(in) :: update_SOC

    logical :: is_KS_basis
    integer(i32) :: kf

    this%IPA = is_IPA
    this%lapwlo_basis = is_LAPWLO_basis
    is_KS_basis = .not. this%lapwlo_basis
    if ( is_KS_basis ) then
      CALL_ASSERT ( max_dimension == n_states, 'max_dimension must be equal to n_states' )
    else
      CALL_ASSERT( max_dimension >= maxval( dims ), 'max_dimension must be >= maxval( dims )' )
    end if
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
    this%obtainSOC = update_SOC
    if ( scissor_shift > eps_scissor .and. is_LAPWLO_basis ) &
      call this%scissor_matrix%allocate_array( [1, 1, ki], [max_dimension, max_dimension, kf] )
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
  subroutine hamiltonian_set_calculate( this, l_max_pot, apwalm, Gkset, &
    ks_lapwlo_transition_matrix, psi_gnd_second_variation, printTimings, t_ham, &
    a_tot, obtain_mathcalH, t_MD ) ! ks_lapwlo_transition_matrix,&
    class(hamiltonian_set), intent(inout) :: this
    !> Maximal value of l in spherical harmonics expansion of DFT potential
    integer(i32), intent(in) :: l_max_pot
    !> Matching coefficients of the (L)APWs (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, :)
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> KS-LAPW+lo transition matrix (nmatmax, n_basis_ks, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: ks_lapwlo_transition_matrix(:, :, :)
    !> Second-variational ground state wavefunctions
    complex(dp), contiguous, optional, intent(in) :: psi_gnd_second_variation(:, :, :)
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

    if ( this%represented_in_lapwlo() ) then
      call this%calculate_in_lapw_basis( l_max_pot, apwalm, Gkset, &
        printTimings, t_ham, a_tot, obtain_mathcalH, t_MD )
    else
      CALL_ASSERT( present( ks_lapwlo_transition_matrix ), "ks_lapwlo_transition_matrix should be present if KS basis is used" )
      call this%calculate_in_ks_basis( l_max_pot, apwalm, Gkset, &
        ks_lapwlo_transition_matrix, psi_gnd_second_variation, printTimings, t_ham )
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
    if( timings_general ) then
      CALL_ASSERT( present( t_ham ) .or. present( t_MD ),  't_ham or t_MD must be present when general timing is desired' )
    end if
    if( timings_detailed ) then
      CALL_ASSERT( timings_general, 'timings_general must be true if timings_detailed is true')
    end if

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
      if ( allocated( this%scissor_matrix%array ) ) this%H_t%array = this%H_t%array + this%scissor_matrix%array
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
      CALL_ASSERT( present( a_tot ), "a_tot must be present" )
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
      CALL_ASSERT( size(gplusk_cart, 1) == n_cartesian, "gplusk_cart must have n_cartesian elements along 1st dim" )
      CALL_ASSERT( size(gplusk_cart, 2) >= n_pw, "gplusk_cart must at least n_pw elements along 2nd dim" )
      CALL_ASSERT( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "mathcalH_ik_ias must have shape [m, m, n_cartesian]" )
      CALL_ASSERT( all(shape(mathcalH_ik_ias) == [m, m, n_cartesian]), "p_MT_ias_ik must have shape [m, m, n_cartesian]" )
      CALL_ASSERT( n_pw <= m, "n_pw must be <= m" )
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
      ks_lapwlo_transition_matrix, psi_gnd_second_variation, printTimings, t_ham )
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
    !> Second-variational ground state wavefunctions
    complex(dp), contiguous, optional, intent(in) :: psi_gnd_second_variation(:, :, :)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(out) :: t_ham

    integer(i32) :: ik, first_kpt, last_kpt, shift, n_basis_first_variation, is, ia, ias
    integer(i32) :: n_atoms, n_states_second_variation, ngp
    complex (dp), allocatable :: V_KS_first_variation(:, :), mt_contribution(:, :, :)
    complex (dp), allocatable :: V_KS_second_variation(:, :), V_SOC(:, :)
    logical :: timings_general, timings_detailed, second_variation, spin_orbit_coupling
    real(dp) :: ti
    real(dp), allocatable :: V_SOC_radial(:, :)

    first_kpt = lbound( this%H_t%array, 3 )
    last_kpt = ubound( this%H_t%array, 3 )
    n_basis_first_variation = size( ks_lapwlo_transition_matrix, 2 )
    shift = first_kpt - 1
    second_variation = present( psi_gnd_second_variation )
    if( second_variation ) then
      n_states_second_variation = size( psi_gnd_second_variation, 1 )
      allocate( V_KS_second_variation(n_states_second_variation, n_states_second_variation) )
      CALL_ASSERT( all( shape( psi_gnd_second_variation ) == shape( this%H_t%array ) ), 'psi_gnd_second_variation must have the same shape as H_t%array' )
    end if
    spin_orbit_coupling = isspinorb() .and. this%obtainSOC
    if( spin_orbit_coupling ) then
      allocate( V_SOC(n_states_second_variation, n_states_second_variation) )
      n_atoms = size( veffmt, 3 )
      allocate( V_SOC_radial(nrcmtmax, n_atoms), source = real_zero )
    end if

    timings_general = .False.
    timings_detailed = .False.
    if ( present( printTimings ) ) call printTimings%get( timings_general, timings_detailed )
    if( timings_general ) then
      CALL_ASSERT( present( t_ham ),  't_ham must be present when general timing is desired' )
    end if
    if( timings_detailed ) then
      CALL_ASSERT( timings_general, 'timings_general must be true if timings_detailed is true')
    end if

    if( timings_general ) call timesec( ti )

    this%H_t%array = zzero
    call this%H_t%copy_from( this%initial_eigenvalues )

    if( this%explicit_evaluation_needed ) then
      call this%H_t%subtract( this%V_KS_0 )

      call me_mt_alloc( mt_contribution )
      allocate( V_KS_first_variation(n_basis_first_variation, n_basis_first_variation) )
    
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! computes gaunts times radial integrals
          call me_mt_prepare( is, ias, lmaxvr, zone, veffmt(:, :, ias), zzero, &
            mt_contribution(:, :, ias) )
        end do ! natoms
      end do ! nspecies
      if( spin_orbit_coupling ) call obtain_SOC_potential_radial( V_SOC_radial )

      do ik = first_kpt, last_kpt
        ngp = Gkset%ngk(1, ik)
        V_KS_first_variation = zzero
        ! mt contribution
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            call me_mt_mat( is, ias, ngp, apwalm(:, :, :, ias, ik-shift), &
              ks_lapwlo_transition_matrix(:, :, ik-shift), zone, mt_contribution(:, :, ias), &
              zone, V_KS_first_variation )
          end do ! natoms
        end do ! nspecies
        ! ir contribution
        call me_ir_mat( Gkset, ik, ks_lapwlo_transition_matrix(:, :, ik-shift), &
          zone, veffig, zone, V_KS_first_variation )
        ! final hamiltonian
        if( second_variation ) then
          call variation_multiplication( psi_gnd_second_variation(:, :, ik-shift), V_KS_first_variation, &
            psi_gnd_second_variation(:, :, ik-shift), V_KS_second_variation, &
            dimA=n_states_second_variation, dimB=n_states_second_variation, startA=1, startB=1 )
          if( spin_orbit_coupling ) then
            call obtain_SOC_potential( lmaxvr, ngp, apwalm(:, :, :, :, ik-shift), &
              ks_lapwlo_transition_matrix(:, :, ik-shift), psi_gnd_second_variation(:, :, ik-shift), V_SOC_radial, V_SOC )
            V_KS_second_variation = V_KS_second_variation + V_SOC
          end if
          this%H_t%array(:, :, ik) = this%H_t%array(:, :, ik) + V_KS_second_variation
        else
          this%H_t%array(:, :, ik) = this%H_t%array(:, :, ik) + V_KS_first_variation
        end if
      end do ! ik
    end if ! this%explicit_evaluation_needed

    if( timings_general ) call timesec_RTTDDFT( ti, t_ham%total )
  end subroutine

  !> (private) Obtain the SOC potential inside each MT sphere \(J\) as
  !> \[ V_{\mathrm{SOC}, J}(r, t) = \frac{\alpha^2/4}{(1-\alpha^2 V_{\mathrm{KS}, J}(r, t)/2)^2} 
  !> \frac{1}{r} \frac{\partial V_{\mathrm{KS}, J}(r,t)}{\partial r}, \]
  !> where \( V_{\mathrm{KS}, J}(r,t) \) is the spherically averaged KS potential inside 
  !> a given MT sphere \(J\), and \( \alpha \) is the fine structure constant. 
  !> If the KS potential inside a MT sphere \(J\) is expanded as
  !> \[ V_{\mathrm{KS}, J}(\mathbf{r}, t) = \sum_{lm} 
  !>          v_{lm, J}(r, t) Y_{lm}(\hat{r}), \]
  !> then, \(V_{\mathrm{KS}, J}(r,t)\) is simply
  !> \[ V_{\mathrm{KS}, J}(r, t) = v_{00,J}(r, t) Y_{00}, \]
  subroutine obtain_SOC_potential_radial( V_SOC_radial )
    !> SOC potential \( V_{\mathrm{KS}, J}(r,t)\) inside each MT sphere. 
    !> 2nd index: MT sphere around each atom; 
    !> 1st index: radial grid of the corresponding MT sphere
    real(dp), contiguous, intent(out) :: V_SOC_radial(:, :)

    integer(i32) :: ia, ias, is, ir, irc, n
    real(dp) :: aux
    real(dp), parameter :: factor1 = 0.5_dp * alpha**2, factor2 = 0.25_dp * alpha**2
    real(dp), allocatable :: vr(:), dv_dr(:), fake(:, :)

    n = size( veffmt, 2 )
    allocate( vr(n), dv_dr(n), fake(3, n) )
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        vr = veffmt(1, :, ias) * y00
        call fderiv(1, nrmt(is), spr(:, is), vr, dv_dr, fake)
        irc = 0
        do ir = 1, nrmt(is), input%groundstate%lradstep
          irc = irc + 1
          aux = 1 - factor1 * vr(ir)
          V_SOC_radial(irc, ias) = factor2 * dv_dr(ir) / (spr(ir, is)*aux**2)
        end do
      end do
    end do
  end subroutine

  !> (private) Obtain the SOC potential in the 2nd variational KS basis at \(t = 0\). 
  !> Inspired in the `seceqnsv` subroutine.
  !> The 2nd variational wavefunctions are assumed to be expanded as
  !> \[ \psi_{i\mathbf{k}}(\mathbf{r}) = \sum_{\sigma n}C^{\mathrm{SV}}_{ni\mathbf{k}\sigma}
  !> \phi_{n\mathbf{k}\sigma}(\mathbf{r})|\sigma\rangle,  \]
  !> where \(\phi_{n\mathbf{k}\sigma}(\mathbf{r})\) are the first variational wavefunctions at \(t = 0\),
  !> \(C^{\mathrm{SV}}_{ni\mathbf{k}\sigma}\) are the second variational coefficients at time \(t = 0\),
  !> and \(|\sigma\rangle\) is the spin state, that could be \(|\uparrow\rangle\) or \(|\downarrow\rangle\). 
  !> The SOC potential determined here can be formally expressed as
  !> \[ \hat{V}_{\mathrm{SOC}}(t) = \sum_{J} V_{\mathrm{SOC}, J}(r, t) \boldsymbol{\sigma}\cdot\mathbf{L}, \]
  !> where \(V_{\mathrm{SOC}, J}(r, t)\) is radial part of the SOC potential inside the MT sphere \(J\)
  !> (as given by [[obtain_SOC_potential_radial]]), 
  !> \(\boldsymbol{\sigma} = \sigma_x \hat{x} + \sigma_y \hat{y} + \sigma_z \hat{z}\) 
  !> is the Pauli matrix, and \(\mathbf{L}\) is the angular momentum operator.
  !> Using the 2nd variational KS wavefunctions at \(t = 0\) as basis:
  !> \[ \hat{V}_{\mathrm{SOC},ij}(t) = \sum_{mn\sigma\sigma'} 
  !> (C^{\mathrm{SV}}_{mi\mathbf{k}\sigma})^* \hat{V}_{\mathrm{SOC},mn}(t)
  !> (C^{\mathrm{SV}}_{nj\mathbf{k}\sigma'}), \]
  !> where
  !> \[ \hat{V}_{\mathrm{SOC},mn}(t) = \sum_J 
  !>    \langle \phi_{m\mathbf{k}\sigma} \sigma| 
  !>     V_{\mathrm{SOC}, J}(r, t) \boldsymbol{\sigma}\cdot\mathbf{L}|
  !>     | \phi_{n\mathbf{k}\sigma'}\sigma\rangle. \] 
  !> This subroutine initially implements the following steps:
  !> <ol>
  !> <li> Given \(J\), the radial components of the 1st variational wavefunctions are explicitly obtained 
  !> \[ \phi_{m\mathbf{k}\sigma}(\mathbf{r})|\sigma\rangle = 
  !>     \sum_{lm} A_{{\bf G+k},lm,\xi} u_{lm}(r) Y_{lm}(\hat{r}) \]
  !> </li> 
  !> <li> The \(\mathbf{L}\) operator is applied to \(Y_{lm}\)</li>
  !> <li> The result is multiplied by \(V_{\mathrm{SOC}, J}(r, t)\) </li>
  !> <li> The \(\boldsymbol{\sigma}\) operator is applied to \(|\sigma'\rangle\) </li>
  !> <li> The result is then used to compute the expectation value with \(\langle \phi_{m\mathbf{k}\sigma} \sigma|\) </li>
  !> <li> The result is accumulated, as a sum over \(m\), \(n\), and \(J\) is required </li>
  !> <li> The matrix operation with \(C^{\mathrm{SV}}_{n\mathbf{k}\sigma'}\) takes place.<\li>
  !> </ol>
  subroutine obtain_SOC_potential( lmax_vr, ngp, apwalm, evecfv, evecsv, V_SOC_radial, V_SOC )
    integer(i32), intent(in) :: lmax_vr
    !> number of \({\bf G+k}\) vectors
    integer(i32), intent(in) :: ngp
    !> wavefunction matching coefficients \(A_{{\bf G+k},lm,\xi}\)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :)
    !> First variation eigenvectors
    complex(dp), contiguous, intent(in) :: evecfv(:, :)
    !> Second variation eigenvectors
    complex(dp), contiguous, intent(in) :: evecsv(:, :)
    !> SOC potential \( V_{\rm SOC}(r)\) inside each MT sphere.
    real(dp), contiguous, intent(in) :: V_SOC_radial(:, :)
    !> SOC potential in the 1st variational KS basis
    complex(dp), contiguous, intent(out) :: V_SOC(:, :)

    integer(i32) :: i, ia, ias, irc, is, ist, j, jst, k, lm, lmmax_vr, n_basis_first_variation, n_r
    real(dp) :: t1
    complex(dp), allocatable :: wfmt1(:, :, :), wfmt2(:, :, :), zlflm(:, :)
    complex(dp), external :: zfmtinp

    V_SOC = zzero
    n_basis_first_variation = size( evecfv, 2 )
    n_r = size( V_SOC_radial, 1 )
    lmmax_vr = ( lmax_vr + 1 ) ** 2
    allocate( wfmt1(lmmax_vr, n_r, n_basis_first_variation), source = zzero )
    allocate( wfmt2(lmmax_vr, n_r, n_cartesian), source = zzero )
    allocate( zlflm(lmmax_vr, n_cartesian), source = zzero )
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        ! radial part of 1st variational wavefunctions is stored in wfmt1
        call generate_basisfunction_secondvariation_MT( lmax_vr, lmmax_vr, ia, is, &
          ngp, apwalm, evecfv, wfmt1 )
        do jst = 1, n_basis_first_variation
          do irc = 1, n_r
            ! apply the \(\mathbf{L}\) operator is applied to each \(lm\)-component
            call lopzflm( lmax_vr, wfmt1(:, irc, jst), lmmax_vr, zlflm )
            t1 = V_SOC_radial(irc, ias)
            ! apply the spin operator
            do lm = 1, lmmax_vr
              ! result for \(\sigma_zL_z|\uparrow\rangle\) is stored in the 1st component
              wfmt2(lm, irc, 1) = wfmt2(lm, irc, 1) + t1 * zlflm(lm, 3)
              ! result for \(\sigma_zL_z|\downarrow\rangle\) is stored in the 2nd component
              wfmt2(lm, irc, 2) = wfmt2(lm, irc, 2) - t1 * zlflm(lm, 3)
              ! result for \((\sigma_xL_x+\sigma_yL_y)|\downarrow\rangle\) is stored in the 3rd component
              wfmt2(lm, irc, 3) = wfmt2(lm, irc, 3) + t1 * ( zlflm(lm, 1) - zi*zlflm(lm, 2) )
            end do
          end do
          do ist = 1, n_basis_first_variation
            do k = 1, n_cartesian
              if (k == 1) then
                ! case: \(\langle \uparrow | \) and \( \uparrow \rangle\)
                i = ist; j = jst
              else if (k == 2) then
                ! case: \(\langle \downarrow | \) and \( \downarrow \rangle\)
                i = ist + n_basis_first_variation; j = jst + n_basis_first_variation
              else
                ! case: \(\langle \uparrow | \) and \( \downarrow \rangle\)
                i = ist; j = jst + n_basis_first_variation
              end if
              V_SOC(i, j) = V_SOC(i, j) + zfmtinp(.True., lmax_vr, n_r, &
                  rcmt(:, is), lmmax_vr, wfmt1(:, :, ist), wfmt2(:, :, k) )
            end do
          end do
        end do
      end do
    end do
    ! Hermitize
    do ist = 1, size( V_SOC, 1 )
      do jst = 1, ist-1
        V_SOC(ist, jst) = conjg( V_SOC(jst, ist) )
      end do
      V_SOC(ist, ist) = V_SOC(ist, ist)%re
    end do
    deallocate( zlflm )
    allocate( zlflm, mold = V_SOC)
    call matrix_multiply( V_SOC, evecsv, zlflm )
    call matrix_multiply( evecsv, zlflm, V_SOC, trans_A='C' )
    CALL_ASSERT( is_hermitian(V_SOC), 'V_SOC is not hermitian' )
  end subroutine
  
  !> Add the pre-calculated length gauge interaction term to the Hamiltonian
  subroutine add_external_coupling_berry_phase( this, external_coupling_length_gauge )
    class(hamiltonian_set), intent(inout) :: this
    !> Length gauge interaction matrix (n_basis, n_basis, n_kpts)
    complex(dp), contiguous, intent(in) :: external_coupling_length_gauge(:, :, :)

    CALL_ASSERT( all( shape( external_coupling_length_gauge ) == shape( this%H_t%array ) ),  'external_coupling_length_gauge and hamiltonian have incompatible dimensions' )

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
    !> Momentum matrix elements
    class(pmat_set), intent(in) :: pmat

    real(dp), parameter :: interaction_tol = 1.e-14_dp
    integer(i32) :: ik, i
    real(dp) :: a_scaled(n_cartesian), fact

    CALL_ASSERT( this%represented_in_lapwlo() .eqv. pmat%represented_in_lapwlo(), "different basis sets used for pmat and Hamiltonian" )
    a_scaled = a_tot%components / c
    fact = 0.5_dp * dot_product( a_scaled, a_scaled )
    if ( fact < interaction_tol ) return
    associate( m => size( this%H_t%array, 1 ), n_kpts => size( this%H_t%array, 3 ) )
      do i = 1, n_cartesian
        CALL_ASSERT( size( pmat%components(i)%array, 3 ) == n_kpts, "pmat(" // to_char( i ) // ") and hamiltonian have different n_kpts" )
      end do
      if( overlap%is_identity() ) then
        ! TODO: use DO CONCURRENT here after ifort2021 support is dropped
        do ik = lbound(this%H_t%array, 3), ubound(this%H_t%array, 3)
          do i = 1, m
            this%H_t%array(i, i, ik) = this%H_t%array(i, i, ik) + cmplx( fact, kind = dp )
          end do
        end do
      else
        CALL_ASSERT( all( shape( overlap%array ) == shape( this%H_t%array ) ),  "overlap and hamiltonian must have same shape" )
        call scaled_add( fact, overlap%array, this%H_t%array )
      end if
      do i = 1, n_cartesian
        call scaled_add( a_scaled(i), pmat%components(i)%array, this%H_t%array )
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

  !> Return whether the Hamiltonian is represented in the LAPW+lo basis
  pure logical function hamiltonian_set_represented_in_lapwlo( this ) result( represented_in_lapwlo )
    class(hamiltonian_set), intent(in) :: this

    represented_in_lapwlo = this%lapwlo_basis
  end function

  !> When the LAPWlo basis is used, build the matrix of the scissor operator \( V_{\rm scissor} \), 
  !> which rigidly shifts the conduction band upwards by \( \Delta E \) to adjust the band gap:
  !> \[
  !> V_{\rm scissor} = \Delta E \sum_{i} {\Theta}(\epsilon_i - E_{\rm Fermi}) \left S | \psi_i \right \rangle 
  !> \left \langle \psi_i \right | S^{\dagger},
  !> \]
  !> where \( \psi_i \) and \( \epsilon_i \) are the ground-state Kohn Sham eigenstates and eigenenergies, 
  !> \( E_{\rm Fermi} \) is the Fermi energy, \( S \) is the overlap matrix, and \( \Theta (x) \) is the Heaviside step function.
  subroutine hamiltonian_set_build_lapwlo_scissor_matrix( this, scissor_shift, &
    first_unoccupied, overlap, ks_lapwlo_transition_matrix )
    class(hamiltonian_set), intent(inout) :: this
    !> Energy shift  \( \Delta E \) for scissor operator
    real(dp), intent(in) :: scissor_shift
    !> Position of the first unoccupied state
    integer, intent(in) :: first_unoccupied
    !> Object that encapsulates the overlap matrix \( S \)
    class(overlap_set), intent(in) :: overlap
    !> KS-LAPW+lo transition matrix (nmatmax, n_ks_states, n_kpt_this_proc)
    complex(dp), contiguous, intent(in):: ks_lapwlo_transition_matrix(:, :, :)

    integer(i32) :: ik, first_kpt, last_kpt, k_shift, n_conduction_states
    complex(dp), allocatable :: overlap_times_psi_lapwlo(:, :)

    if ( scissor_shift > eps_scissor ) then
      call this%scissor_matrix%assert_allocated()
      CALL_ASSERT ( size( this%initial_eigenvalues, 1 ) == size( ks_lapwlo_transition_matrix, 2 ), '2nd dimension of ks_lapwlo_transition_matrix must be equal to the 1st dimension of initial_eigenvalues' )
      CALL_ASSERT ( size( this%H_t%array, 3 ) == size( ks_lapwlo_transition_matrix, 3 ), '3rd dimension of ks_lapwlo_transition_matrix and this%H_t%array must be the same' )
      n_conduction_states = size( this%initial_eigenvalues, 1 ) - first_unoccupied + 1
      first_kpt = lbound( this%H_t%array, 3 )
      last_kpt = ubound( this%H_t%array, 3 )
      allocate( overlap_times_psi_lapwlo(size( ks_lapwlo_transition_matrix, 1 ), n_conduction_states), source = zzero )
      k_shift = 1 - first_kpt
      do ik = first_kpt, last_kpt
        overlap_times_psi_lapwlo = zzero
        call matrix_multiply( overlap%array(:, :, ik ), &
          ks_lapwlo_transition_matrix(:, first_unoccupied:, ik + k_shift), overlap_times_psi_lapwlo )
        call matrix_multiply( overlap_times_psi_lapwlo, overlap_times_psi_lapwlo, &
          this%scissor_matrix%array(:, :, ik), trans_B = 'C' )
        this%scissor_matrix%array(:, :, ik) = this%scissor_matrix%array(:, :, ik) * scissor_shift
      end do
    end if
  end subroutine

end module rttddft_Hamiltonian
