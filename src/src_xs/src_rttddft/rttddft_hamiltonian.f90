! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Reference: https://doi.org/10.1088/2516-1075/ac0c26
! TODO(Ronaldo): Refactor to reduce the number of global variables
!> Module that manages the hamiltonian matrix in RT-TDDFT
module rttddft_Hamiltonian
  use asserts, only: assert
  use constants, only: fourpi, y00, zi, zone, zzero
  use matrix_elements, only: me_mt_alloc, me_mt_prepare, me_mt_mat, me_ir_alloc, me_ir_prepare, me_ir_mat
  use mod_atoms, only: atposc, idxas, natoms, nspecies
  use mod_gvector, only: cfunig
  use mod_kpointset, only: Gk_set
  use mod_lattice, only: omega
  use mod_muffin_tin, only: rmt
  use mod_potential_and_density, only: meffig, veffig, veffir, veffmt
  use modinput, only: input
  use physical_constants, only: alpha, c
  use precision, only: dp, i32
  use rttddft_GlobalMDVariables, only: mathcalH
  use rttddft_Overlap, only: overlap_set
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_hamiltonian, &
    Timing_Ehrenfest, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential_Field

  implicit none

  private
  public :: update_hamiltonian_without_pa_term_lapw, &
    add_external_coupling_berry_phase, update_hamiltonian_without_pa_term_ks, &
    add_external_coupling_velocity_gauge
  
  integer(i32), parameter :: n_cartesian = 3

contains

  !> In `update_hamiltonian_without_pa_term_lapw`, we obtain the explicitly 
  !> field-free hamiltonian at time \( t \) in the LAPW+lo basis.
  subroutine update_hamiltonian_without_pa_term_lapw( first_kpt, l_max_pot, a_tot, ham_time, apwalm, Gkset, &
    printTimings, t_ham, t_MD, update_mathcalH )
    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Maximal value of l in spherical harmonics expansion of DFT potential
    integer(i32), intent(in) :: l_max_pot
    !> Total vector potential
    type(Vector_Potential_Field), intent(in) :: a_tot
    !> Hamiltonian matrix at current time \(t\) (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), contiguous, intent(inout) :: ham_time(:, :, first_kpt :)
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, first_kpt :)
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(out) :: t_ham
    !> Object that packs information about timings related to MD
    type(Timing_Ehrenfest), optional, intent(out) :: t_MD
    !> if `.True.`, update `mathcalH`
    logical, intent(in), optional :: update_mathcalH

    integer(i32) :: ik, last_kpt, ias, ia, is, n_MT_radial_points
    real(dp) :: t_i, t_f, t_aux
    logical :: timings_general, timings_detailed, get_mathcalH, valence_relativity
    integer(i32), parameter :: l_max_kin = 0 ! Overlap operator is equal to (1.0/y00)*Y_{00}, it has only l=0 component
    integer(i32), parameter :: lm_max = (l_max_kin+1)**2
    real(dp), parameter :: aux_factor = 0.5_dp*alpha**2*y00
    real(dp), parameter :: kin_operator_00 = (0.5_dp / y00)
    real(dp), allocatable :: kin_mt(:, :)
    complex(dp), allocatable :: mt_part(:, :, :)

    last_kpt = ubound( ham_time, 3 )
    ! Check optional arguments
    timings_general = .False.
    timings_detailed = .False.
    if ( present( printTimings ) ) call printTimings%get( timings_general, timings_detailed )
    get_mathcalH = .False.
    if( present( update_mathcalH ) ) get_mathcalH = update_mathcalH

    ! sanity checks
    if( timings_general ) call assert( present( t_ham ) .or. present( t_MD ), &
      't_ham or t_MD must be present when general timing is desired' )
    if( timings_detailed ) call assert( timings_general, 'timings_general must be true if timings_detailed is true')

    if( timings_general ) call timesec( t_i )

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
        call calculate_hamiltonian_without_pa_term_lapw_ik( ik, Gkset, apwalm(:, :, :, :, ik), &
          mt_part, meffig, ham_time(:, :, ik), get_mathcalH, a_tot%components )
      end do
    else
      do ik = first_kpt, last_kpt
        call calculate_hamiltonian_without_pa_term_lapw_ik( ik, Gkset, apwalm(:, :, :, :, ik), &
          mt_part, cfunig, ham_time(:, :, ik), get_mathcalH, a_tot%components )
      end do
    end if

    if( timings_general ) then
      call timesec( t_f )
      if( present( t_ham ) ) t_ham%total = t_f - t_i
      if( timings_detailed .and. present( t_MD ) ) t_MD%ham = t_f - t_i
    end if
  end subroutine update_hamiltonian_without_pa_term_lapw

  !> Calculate the hamiltonian matrix for a given k-point.
  subroutine calculate_hamiltonian_without_pa_term_lapw_ik( ik, Gkset, apwalm_ik, mt_part, &
      kin_ir, ham_ik, get_mathcalH, a_tot )
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
    !> MT part of the overlap matrix
    complex(dp), contiguous, intent(inout) :: ham_ik(:, :)
    !> If `.True.`, calculate the MT contributions to the auxiliary matrix mathcalH
    logical, intent(in) :: get_mathcalH
    !> `x`, `y`, and `z` components of the (total) vector potential
    real(dp), optional, intent(in) :: a_tot(n_cartesian)

    integer(i32) :: ias, ia, is
    complex(dp), allocatable :: tmp(:, :)

    ham_ik = zzero
    ! It is better to split the two cases (even with some code duplication): 
    ! - to avoid an "if" inside the double loop over atoms
    ! - to avoid allocating "tmp" when not needed (this can be a large array)
    associate( np => Gkset%ngk(1, ik), m => size(ham_ik, 1) )
    if( get_mathcalH ) then
      tmp = ham_ik
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, np, apwalm_ik(:, :, :, ias), zone, &
            mt_part(:, :, ias), zzero, tmp )
          ham_ik = ham_ik + tmp
          call update_mathcalH_ik_ias( is, ia, np, Gkset%vgkc(:, :, 1, ik), tmp, a_tot, mathcalH(:, :, :, ias, ik) )
        end do
      end do
    else
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, np, apwalm_ik(:, :, :, ias), zone, &
            mt_part(:, :, ias), zone, ham_ik )
        end do
      end do
    end if
    end associate
    call me_ir_mat( Gkset, ik, zone, veffig, zone, ham_ik )
    call me_ir_mat( Gkset, ik, zone/2, kin_ir, zone, ham_ik, gradient_product=.true. )
  end subroutine

  !> (private) Update the `mathcalH` matrix for a given atom and k-point
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

  !> In `update_hamiltonian_without_pa_term_ks`, we obtain the explicitly field-independent 
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
  subroutine update_hamiltonian_without_pa_term_ks( first_kpt, lmaxvr, ham_time, apwalm, &
      ks_lapwlo_transition_matrix, effective_potential_init, ham_init, Gkset, printTimings, t_ham )
    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Maximal value of l in spherical harmonics expansion of DFT potential
    integer(i32), intent(in) :: lmaxvr
    !> Field-free Hamiltonian matrix \( H(t) \) (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), contiguous, intent(out) :: ham_time(:, :, first_kpt :)
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, first_kpt :)
    !> KS-LAPW+lo transition matrix (nmatmax, n_basis_ks, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: ks_lapwlo_transition_matrix(:, :, first_kpt :)
    !> Initial effective potential \( V_{\rm eff}(t = 0) \) (n_basis, n_basis, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: effective_potential_init(:, :, first_kpt :)
    !> Hamiltonian matrix \( H_{\rm init} \) obtained at time \(t = 0 \)
    complex(dp), contiguous, intent(in) :: ham_init(:, :, first_kpt :)
    !> Set of G+k vectors for LAPW expansion
    type(Gk_set), intent(in) :: Gkset
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(out) :: t_ham

    integer(i32) :: ik, last_kpt, n_basis, is, ia, ias, ngp
    complex (dp), allocatable :: local_effective_potential(:, :), mt_contribution(:, :, :)
    logical :: timings_general, timings_detailed
    real(dp) :: ti

    last_kpt = ubound( ham_time, 3 )
    n_basis = size( ham_time, 1 )

    timings_general = .False.
    timings_detailed = .False.
    if ( present( printTimings ) ) call printTimings%get( timings_general, timings_detailed )
    if( timings_general ) call assert( present( t_ham ), &
      't_ham must be present when general timing is desired' )
    if( timings_detailed ) call assert( timings_general, 'timings_general must be true if timings_detailed is true')

    if( timings_general ) call timesec( ti )

    call me_mt_alloc( mt_contribution )
    allocate( local_effective_potential(n_basis, n_basis), source = zzero )
    
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
               
          call me_mt_mat( is, ias, ngp, apwalm(:, :, :, ias, ik), &
            ks_lapwlo_transition_matrix(:, :, ik), zone, mt_contribution(:, :, ias), &
            zone, local_effective_potential )

        end do ! natoms
      end do ! nspecies

      ! ir contribution
      call me_ir_mat( Gkset, ik, Gkset, ik, ks_lapwlo_transition_matrix(:, :, ik), &
        ks_lapwlo_transition_matrix(:, :, ik), zone, veffig, zone, local_effective_potential )

      ham_time(:, :, ik) = ham_init(:, :, ik) + local_effective_potential - effective_potential_init(:, :, ik)

    end do ! ik

    if( timings_general ) call timesec_RTTDDFT( ti, t_ham%total )

  end subroutine
  
  !> Add the pre-calculated length gauge interaction term to the Hamiltonian
  subroutine add_external_coupling_berry_phase( external_coupling_length_gauge, ham_time, dims )
    !> Length gauge interaction matrix (n_basis, n_basis, n_kpts)
    complex(dp), contiguous, intent(in) :: external_coupling_length_gauge(:, :, :)
    !> Hamiltonian matrix at current time \(t\) (n_basis, n_basis, n_kpts)
    complex(dp), contiguous, intent(inout) :: ham_time(:, :, :)
    !> Used dimensions of `external_coupling_length_gauge` and `ham_time` matrices
    integer(i32), intent(in) :: dims(:)

    integer(i32) :: ik

    call assert( all( shape( external_coupling_length_gauge ) == shape( ham_time ) ), &
      'external_coupling_length_gauge and ham_time have incompatible dimensions' )

    !$omp parallel default(none), private(ik), &
    !$omp& shared(ham_time, external_coupling_length_gauge, dims)
    !$omp do
    do ik = 1, size( ham_time, 3 )
      ham_time(1 : dims(ik), 1 : dims(ik), ik) = ham_time(1 : dims(ik), 1 : dims(ik), ik) + &
        external_coupling_length_gauge(1 : dims(ik), 1 : dims(ik), ik)
    end do
    !$omp end parallel

  end subroutine

  !> Add the velocity gauge interaction term \( {\bf p} \cdot {\bf A}(t) / c \) to 
  !> the Hamiltonian at time \( t \).
  ! TODO: is the space-uniform A^2 term needed here?
  subroutine add_external_coupling_velocity_gauge( a_tot, overlap, ham_time, pmat, dims )
    !> Total vector potential
    type(Vector_Potential_Field), intent(in) :: a_tot
    !> Overlap matrix
    class(overlap_set), intent(in) :: overlap
    !> Hamiltonian matrix at current time \(t\) (n_basis, n_basis, n_kpts_kpt)
    complex(dp), contiguous, intent(inout) :: ham_time(:, :, :)
    !> Momentum matrix elements (n_basis, n_basis, 3, n_kpts)
    complex(dp), contiguous, intent(in) :: pmat(:, :, :, :)
    !> Used dimensions of `overlap` and `ham_time` matrices
    integer(i32), intent(in) :: dims(:)

    real(dp), parameter :: interaction_tol = 1.e-14_dp
    integer(i32) :: ik, i
    real(dp) :: a_scaled(3), fact

    a_scaled = a_tot%components / c
    fact = 0.5_dp * dot_product( a_scaled, a_scaled )
    if ( fact < interaction_tol ) return
    associate( m => size( ham_time, 1 ), n_kpts => size( ham_time, 3 ) )
      call assert( size( dims, 1 ) == n_kpts, "dims and ham_time have different n_kpts" )
      call assert( size( pmat, 4 ) == n_kpts, "pmat and ham_time have different n_kpts" )
      if( overlap%is_identity() ) then
        do concurrent( i = 1 : m, ik = 1 : n_kpts )
          ham_time(i, i, ik) = ham_time(i, i, ik) + fact
        end do
      else
        call assert( all( shape( overlap%array ) == shape( ham_time ) ), &
          "overlap and ham_time must have same shape" )
        ! TODO: replace by zaxpy wrapper after MR 840 is merged
        ham_time = ham_time + fact * overlap%array
      end if
      ! TODO: replace by zaxpy wrapper after MR 840 is merged
      ham_time = ham_time + a_scaled(1)*pmat(:, :, 1, :) + &
                          + a_scaled(2)*pmat(:, :, 2, :) + &
                          + a_scaled(3)*pmat(:, :, 3, :)
    end associate
  end subroutine

end module rttddft_Hamiltonian
