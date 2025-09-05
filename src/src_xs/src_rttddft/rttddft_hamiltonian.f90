! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Reference: https://doi.org/10.1088/2516-1075/ac0c26
! TODO(Ronaldo): Refactor to reduce the number of global variables
!> Module that manages the hamiltonian matrix in RT-TDDFT
module rttddft_Hamiltonian
  use asserts, only: assert
  use constants, only: fourpi, zi, zone, zzero
  use mod_APW_LO, only: apword, nlorb, lorbl
  use mod_atoms, only: nspecies, natoms, idxas, atposc
  use mod_eigensystem, only: nmat, idxlo, h1aa, h1loa, h1lolo, &
    oalo, ololo, MTHamiltonianList, MTInitAll, MTNullify, MaxAPWs
  use mod_gkvector, only: ngk, vgkc, igkig
  use mod_gvector, only: ivg, ivgig, cfunig, ngvec
  use mod_lattice, only: omega
  use mod_muffin_tin, only: idxlm, rmt
  use mod_potential_and_density, only: veffig, meffig, m2effig, veffmt
  use modinput, only: input
  use physical_constants, only: alpha, c
  use precision, only: dp, i32
  use rttddft_GlobalMDVariables, only: mathcalH, mathcalB
  use rttddft_Overlap, only: overlap_set
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_hamiltonian, &
    Timing_Ehrenfest, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use mod_kpointset, only: Gk_set
  use matrix_elements, only: me_mt_alloc, me_mt_prepare, me_mt_mat, me_ir_mat
  use xlapack, only: matrix_multiply

  implicit none

  private
  public :: update_hamiltonian_without_pa_term_lapw, &
    add_external_coupling_berry_phase, update_hamiltonian_without_pa_term_ks, &
    add_external_coupling_velocity_gauge
    
  real(dp) :: fact, atot(3)
  type(MTHamiltonianList) :: mt_h

contains

  !> In `update_hamiltonian_without_pa_term_lapw`, we obtain the explicitly 
  !> field-free hamiltonian at time \( t \) in the LAPW+lo basis.
  subroutine update_hamiltonian_without_pa_term_lapw( first_kpt, a_tot, ham_time, apwalm, &
    printTimings, t_ham, t_MD, update_mathcalH )
    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Total vector potential
    type(Vector_Potential_Field), intent(in) :: a_tot
    !> Hamiltonian matrix at current time \(t\) (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), contiguous, intent(inout) :: ham_time(:, :, first_kpt :)
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: apwalm(:, :, :, :, first_kpt :)
    !> Object that packs information about printing of timings
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings to update the Hamiltonian
    type(Timing_RTTDDFT_hamiltonian), optional, intent(out) :: t_ham
    !> Object that packs information about timings related to MD
    type(Timing_Ehrenfest), optional, intent(out) :: t_MD
    !> if `.True.`, update `mathcalH`
    logical, intent(in), optional :: update_mathcalH

    integer(i32) :: ik, last_kpt
    real(dp) :: ti, tf, tStart
    logical :: timings_general, timings_detailed, get_mathcalH

    last_kpt = ubound( ham_time, 3 )
  
    ! atot and fact are used in `obtain_interstitial_contribution_mathcalH` through the module
    atot = a_tot%components
    fact = dot_product( atot, atot ) / (2._dp * c**2)

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

    if( timings_general ) then 
      call timesec( ti )
      tStart = ti
    end if

    call mt_h%release()
    call MTNullify( mt_h )
    call MTInitAll( mt_h )
    call hmlint( mt_h )
    if ( timings_detailed .and. present( t_ham ) ) call timesec_RTTDDFT( ti, t_ham%hmlint )

    !$omp parallel default(none), private(ik), &
    !$omp& shared(first_kpt, last_kpt, ham_time, apwalm, nmat, get_mathcalH)
    !$omp do
    do ik = first_kpt, last_kpt
      call hamsetup( ik, ham_time(:, :, ik), apwalm(:, :, :, :, ik), nmat(1, ik), get_mathcalH )
    end do
    !$omp end do 
    !$omp end parallel

    if ( get_mathcalH ) call obtain_interstitial_contribution_mathcalH( &
      & first_kpt, last_kpt )

    if( timings_general ) then
      call timesec( tf )
      if( present( t_ham ) ) t_ham%total = t_ham%total + tf - tStart
      if( timings_detailed .and. present( t_ham ) ) t_ham%rest = tf - ti
      if( timings_detailed .and. present( t_MD ) ) t_MD%ham = t_MD%ham + tf - ti
    end if

  end subroutine update_hamiltonian_without_pa_term_lapw

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

    integer(i32) :: ik, nmatp, last_kpt, apwordmax, lmmaxapw, &
      n_basis, is, ia, ias, ngp
    complex (dp), allocatable :: local_effective_potential(:, :), mt_contribution(:, :, :)
    logical :: timings_general, timings_detailed
    real(dp) :: ti

    last_kpt = ubound( ham_time, 3 )
    apwordmax = size( apwalm, 2 )
    lmmaxapw = size( apwalm, 3 )
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
      ngp = ngk(1, ik)
      nmatp = nmat(1, ik)
      
      local_effective_potential = zzero

      ! mt contribution
      do is = 1, nspecies
        do ia = 1, natoms(is)

          ias = idxas(ia, is)
               
          call me_mt_mat( is, ias, ngp, ngp, apwalm(1 : ngp, :, :, ias, ik), &
            apwalm(1 : ngp, :, :, ias, ik), ks_lapwlo_transition_matrix(1 : nmatp, :, ik), &
            ks_lapwlo_transition_matrix(1 : nmatp, :, ik), zone, mt_contribution(:, :, ias), &
            zone, local_effective_potential )

        end do ! natoms
      end do ! nspecies

      ! ir contribution
      call me_ir_mat( Gkset, ik, Gkset, ik, ks_lapwlo_transition_matrix(1 : nmatp, :, ik), &
        ks_lapwlo_transition_matrix(1 : nmatp, :, ik), zone, veffig, zone, local_effective_potential )

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
    real(dp) :: a_scaled(3)

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

  !> Subroutine to calculate the interstitial contribution to `mathcalH` (used to
  !> obtain the forces on the ions in Ehrenfest Dynamics)
  subroutine obtain_interstitial_contribution_mathcalH( first_kpt, last_kpt )
    integer, intent(in) :: first_kpt, last_kpt

    integer     :: ik, is, ia, ias, ig, igl, j, ngp
    real(dp)    :: t1, t2, t3, t4, g(3)
    complex(dp) :: t5

    !$omp parallel default(none), private(ik,ngp,is,ia,ias,ig,igl,g,t1,t2,t3,t4,t5), &
    !$omp& shared(first_kpt,last_kpt,natoms), &
    !$omp& shared(nspecies,rmt,omega,atposc,idxas,atot), &
    !$omp& shared(ngk,nmat,vgkc,igkig,mathcalH)
    !$omp do
    do ik = first_kpt, last_kpt
      ngp = ngk(1,ik)
      ! Loop over atoms
      do is = 1, nspecies
        t1 = fourpi*(rmt(is)**3)/omega
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! Loop over g-points
          do ig = 1, ngp
            do igl = 1, ngp
              if ( ig .eq. igl ) cycle
              t2 = dot_product(0.5d0*vgkc(:,igl,1,ik)+atot(:)/c,vgkc(:,ig,1,ik))
              g(:) = vgkc(:,ig,1,ik)-vgkc(:,igl,1,ik)
              t3 = rmt(is)*dsqrt(g(1)**2+g(2)**2+g(3)**2)
              ! Spherical Bessel Function of 1st kind over t3
              t3 = (sin(t3)-t3*cos(t3))/(t3**3)
              t4 = dot_product(g(:),atposc(:, ia, is))
              t5 = cmplx(cos(t4),sin(t4),kind(dp))
              ! Loop over cartesian coordinates
              do j = 1, 3
              mathcalH(igl,ig,j,ias,ik) = mathcalH(igl,ig,j,ias,ik) - &
                & zi*(g(j))*t1*t2*t3*t5
              end do ! do j = 1, 3
            end do ! do ig = 1, ngp
          end do ! do igl = 1, ngp
        end do ! do ia = 1, natoms(is)
      end do ! do is = 1, nspecies
    end do
    !$omp end do 
    !$omp end parallel
  end subroutine obtain_interstitial_contribution_mathcalH

  !> Subroutine to calculate the Hamiltonian matrix in the LAPW+lo basis for a given k-point
  subroutine hamsetup( ik, ham_time, apwalm, nmatp, calculate_mathcalH )

    !> ik: the index of the k-point considered
    integer, intent(in) :: ik
    !> Hamiltonian matrix at current time \(t\) at the current k-point (nmatmax, nmatmax)
    complex(dp), intent(out) :: ham_time(:, :)
    !> Matching coefficients of the (L)APWs at the current k-point
    !> (ngkmax, apwordmax, lmmaxapw, natmtot)
    complex(dp), intent(in) :: apwalm(:, :, :, :)
    !> `nmatp` is the dimension of the matrix for this `k-point` 
    !> (`nmatp` \( \times \) `nmatp`)
    integer, intent(in) :: nmatp
    !> Is it required to calculate the MT contributions to the auxiliary matrix mathcalH
    logical, intent(in) :: calculate_mathcalH
    ! TODO: calling hamsetup with calculate_mathcalH = .true. without 
    ! calling overlapsetup in advance makes no sense, add check

    integer                   :: i, j, is, ia, ias, if3, ig, igl, io2
    integer                   :: j1, l3, m3, lm3, j3, maxnlo, maxaa
    integer                   :: iv(3)
    integer                   :: ngp
    real (dp)                 :: t1
    complex (dp)              :: zt
    complex (dp), allocatable :: hamcopy(:, :), aux(:, :)
    complex (dp), allocatable :: apwi(:, :), zm(:, :), aux_mathcalh(:, :, :)

    ! auxiliary variables
    ngp = ngk(1,ik)
    maxaa = mt_h%maxaa
    maxnlo = mt_h%maxnlo
    allocate( apwi(maxaa, ngp) )
    allocate( aux(nmatp, nmatp) )
    allocate( zm(maxaa, ngp) )
    allocate( hamcopy(nmatp, nmatp) )
    hamcopy(:,:) = zzero
    if ( calculate_mathcalH ) allocate( aux_mathcalh(nmatp, nmatp, 3) )
    do is = 1, nspecies
      do ia = 1, natoms(is)
        if ( calculate_mathcalH ) aux_mathcalh = zzero
        ! APW-APW part
        ias = idxas (ia, is)
        apwi = zzero
        if3 = 0
        do l3 = 0, input%groundstate%lmaxmat
          do m3 = -l3, l3
          lm3 = idxlm (l3, m3)
            do io2 = 1, apword (l3, is)
              if3 = if3 + 1
              apwi(if3,1:ngp) = apwalm(1:ngp,io2,lm3,ias)
            end do
          end do
        end do
        zm(:,:) = zzero
        ! Matrix multiplication: zm = (muffintin_hamiltonian)*(matching coefficients)
        ! zm = (mt_h%maxaa)*(apwi)
        call ZGEMM( 'N', 'N', maxaa, ngp, maxaa, zone, &
          & mt_h%main%aa(:,:,ias), maxaa, apwi, maxaa, zone, zm, maxaa )
        ! Matrix multiplication: hamcopy = hamcopy + (matching coefficients)^H*(zm)
        call ZGEMM( 'C', 'N', ngp, ngp, maxaa, zone, apwi, maxaa, zm, maxaa, &
          & zzero, aux, nmatp )
        hamcopy(1:ngp, 1:ngp) = hamcopy(1:ngp, 1:ngp) + aux(1:ngp, 1:ngp)
        if ( calculate_mathcalH ) then
          do ig = 1, ngp
            do igl = 1, ngp
              aux_mathcalh(igl, ig, 1 : 3) = zi*( vgkc(1:3,ig,1,ik)-vgkc(1:3,igl,1,ik) )*aux(igl,ig)
              mathcalH(igl,ig,1:3,ias,ik) = mathcalH(igl,ig,1:3,ias,ik) + aux_mathcalh(igl, ig, 1 : 3) 
            end do
          end do
        end if

    !What if it is, say, LAPW calculation without any local orbitals?
        if ( nlorb(is) /= 0 ) then
    ! APW-LO part
          l3 = lorbl(1,is)
          lm3 = idxlm(l3,-l3)
          call ZGEMM( 'N', 'N', mt_h%losize(is), ngp, maxaa, zone, &
            & mt_h%main%loa(:,:,ias), maxnlo, apwi, maxaa, zzero, &
            & aux(ngp+idxlo(lm3,1,ias),1), nmatp )
          j1 = ngp + idxlo( lm3, 1, ias )
          j3 = j1 + mt_h%losize(is) - 1
          hamcopy(j1:j3,1:ngp) = hamcopy(j1:j3,1:ngp) + aux(j1:j3,1:ngp)
          if ( calculate_mathcalH ) then
            do ig = 1, ngp
              do j = 1, 3
                aux_mathcalh(j1:j3, ig, j) = zi*(vgkc(j,ig,1,ik))*aux(j1:j3,ig)
                mathcalH(j1:j3,ig,j,ias,ik) = mathcalH(j1:j3,ig,j,ias,ik) + aux_mathcalh(j1:j3, ig, j)                  
              end do
            end do
          end if
          do i = j1, j3
            hamcopy(1:ngp,i)=conjg(hamcopy(i,1:ngp))
            if ( calculate_mathcalH ) then
              do j = 1, 3
                mathcalH(1:ngp,i,j,ias,ik) = mathcalH(1:ngp,i,j,ias,ik) + conjg(aux_mathcalh(i,1:ngp,j))
              end do
            end if
          enddo
    ! LO-LO part
          hamcopy(j1:j3,j1:j3) = hamcopy(j1:j3, j1:j3) + &
            & mt_h%main%lolo(1:1+j3-j1, 1:1+j3-j1, ias)
        endif
      end do
    end do

    ! interstitial contributions
    if ( input%groundstate%ValenceRelativity /= "none" ) then
      do j = 1, ngp
        do i = 1, j
          iv(:) = ivg(:,igkig(i,1,ik)) - ivg(:,igkig(j,1,ik))
          ig = ivgig(iv(1),iv(2),iv(3))
          if ((ig .gt. 0) .and. (ig .le. ngvec)) then
            t1 = 0.5_dp*dot_product(vgkc(:,i,1,ik),vgkc(:,j,1,ik))
            zt = veffig(ig) + t1*meffig(ig)
            hamcopy(i,j) = hamcopy(i,j) + zt
            hamcopy(j,i) = conjg(hamcopy(i,j))
          end if ! if ((ig .gt. 0) .and. (ig .le. ngvec))
        end do ! do i = 1, j
      end do ! do j = 1, ngp
    else
      do j = 1, ngp
        do i = 1, j
          iv(:) = ivg(:,igkig(i,1,ik)) - ivg(:,igkig(j,1,ik))
          ig = ivgig(iv(1),iv(2),iv(3))
          if ((ig .gt. 0) .and. (ig .le. ngvec)) then
            t1 = 0.5_dp*dot_product(vgkc(:,i,1,ik),vgkc(:,j,1,ik))
            zt = veffig(ig) + t1*cfunig(ig)
            hamcopy(i,j) = hamcopy(i,j) + zt
            hamcopy(j,i) = conjg(hamcopy(i,j))
          end if ! if ((ig .gt. 0) .and. (ig .le. ngvec))
        end do ! do i = 1, j
      end do ! do j = 1, ngp
    endif

    ham_time(1 : nmatp, 1 : nmatp) = hamcopy(1 : nmatp, 1 : nmatp)
    
  end subroutine hamsetup

end module rttddft_Hamiltonian
