! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created Apr 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module that manages the hamiltonian and overlap in RT-TDDFT
!> //TODO(Ronaldo): Refactor to reduce the number of global variables
module rttddft_HamiltonianOverlap
  use asserts, only: assert
  use constants, only: fourpi, zi
  
  use mod_APW_LO, only: apword, nlorb, lorbl
  use mod_atoms, only: nspecies, natoms, idxas, natmtot, atposc
  use mod_eigensystem, only: nmat, nmatmax, hloloij, idxlo, h1aa, h1loa, h1lolo, &
    oalo, ololo, MTHamiltonianList, MTInitAll, MTNullify, MTRelease
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_gkvector, only: ngk, vgkc, igkig
  use mod_gvector, only: ivg, ivgig, cfunig, ngvec
  use mod_kpoint, only: nkpt
  use mod_lattice, only: omega
  use mod_muffin_tin, only: idxlm, rmt
  use mod_potential_and_density, only: veffig, meffig, m2effig
  use modinput, only: input
  use modmpi
  use rttddft_GlobalVariables, only: ham_time, ham_past, overlap, mathcalH, apwalm, &
    & atot, pmat, pmatmt, evecfv_time, timesecRTTDDFT
  use rttddft_pmat, only: Obtain_Pmat_LAPWLOBasis
  use physical_constants, only: c
  use precision, only: dp
  
  implicit none

  private
  public :: UpdateHam

  real(dp)              :: fact
  type(MTHamiltonianList) :: mt_h

contains

  !> In UpdateHam, we obtain the hamiltonian (and if requested, the overlap) at 
  !> time \( t \).
  subroutine UpdateHam( predcorr, calculateOverlap, &
    & timeGen, timeDetail, timeini, timefinal, tgenpmatbasis, thmlint, tham, &
    & update_mathcalH, update_mathcalB, update_pmat )
    !> tells if we are in the loop of the predictor-Corrector scheme    
    logical, intent(in)               :: predcorr
    !> tells if we need to calculate the overlap
    logical, intent(in)               :: calculateOverlap
    !> tells if we want a general timing
    logical, intent(in), optional     :: timeGen
    !> tells if we want a detailed timing (only works if `timeGen` is true)
    logical, intent(in), optional     :: timeDetail
    !> time (in seconds) when the subroutine was called (used for making time differences)
    real(dp), intent(in), optional    :: timeini
    !> time (in seconds) after executing this subroutine
    real(dp), intent(out), optional   :: timefinal
    !> time spent to execute hmlint
    real(dp), intent(out), optional   :: thmlint
    !> time spent after executing hmlint until the end of this subroutine
    real(dp), intent(out), optional   :: tham
    !> time spent when generating `pmat`
    real(dp), intent(out), optional   :: tgenpmatbasis
    !> if `.True.`, update `mathcalH`
    logical, intent(in), optional     :: update_mathcalH
    !> if `.True.`, update `mathcalB`
    logical, intent(in), optional     :: update_mathcalB
    !> if `.True.`, update `pmat`
    logical, intent(in), optional     :: update_pmat

    integer               :: ik, nmatp, first_kpt, last_kpt
    real(dp)              :: timei, timef
    logical               :: tGen, tDetail
    logical               :: get_mathcalH, get_mathcalB, get_pmat, forcePmatHermitian

    ! factor that multiplies the overlap matrix (when we compute the hamiltonian)
    fact = dot_product( atot, atot )/(2._dp * c**2)
    

    ! Check optional arguments
    tGen = .False.
    tDetail = .False.
    if ( present(timeGen) ) then
      tGen = timeGen
      if ( present(timeDetail) ) tDetail = timeDetail
    end if
    get_mathcalH = .False.
    if( present( update_mathcalH ) ) get_mathcalH = update_mathcalH
    get_mathcalB = .False.
    if( present( update_mathcalB ) ) get_mathcalB = update_mathcalB
    get_pmat = .False.
    if( present( update_pmat ) ) get_pmat = update_pmat
    forcePmatHermitian = input%xs%realTimeTDDFT%forcePmatHermitian

    ! sanity checks
    if( get_mathcalH ) call assert( calculateOverlap , 'The overlap matrix is needed to update mathcalH' )
    if( tGen ) call assert( present(timeini) .and. present(timefinal), &
      'timeini and timefinal must be present when general timing is desired' )
    if( tDetail ) then 
      call assert( present(thmlint) .and. present(tham), &
      'thmlint tham must be present when detailed timing is desired')
      if( get_pmat ) call assert( present(tgenpmatbasis), &
        'tgenpmatbasis must be present if detailed timing is desired and pmat is updated')
    end if

    if( tGen ) timei = timeini
    if( get_pmat ) then
      call Obtain_Pmat_LAPWLOBasis( forcePmatHermitian, allocated(pmatmt) )
      if( tDetail ) call timesecRTTDDFT( timei, timef, tgenpmatbasis )
    end if

    call MTNullify(mt_h)
    call MTInitAll(mt_h)
    call hmlint(mt_h)

    if ( tDetail ) call timesecRTTDDFT( timei, timef, thmlint )

    if ( .not. predcorr ) ham_past(:,:,:) = ham_time(:,:,:)

    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)
 
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik,nmatp), &
!$OMP& SHARED(first_kpt,last_kpt,calculateOverlap,nkpt,natmtot,natoms), &
!$OMP& SHARED(nspecies,rmt,omega,atposc,idxas,apwalm,fact,atot,pmat,ham_time), &
!$OMP& SHARED(overlap,ngk,nmat,vgkc,igkig,input, get_mathcalH, get_mathcalB), &
!$OMP& SHARED(nmatmax)
!$OMP DO
#endif
    do ik = first_kpt, last_kpt
      nmatp = nmat(1,ik)
      call hamsetup( ik, nmatp, get_mathcalH )
      if ( calculateOverlap ) then
        call overlapsetup( ik, nmatp, get_mathcalB, get_mathcalH )
      end if

      ! Include the part of the vector potential in the hamiltonian
      ham_time(1:nmatp,1:nmatp,ik) = ham_time(1:nmatp,1:nmatp,ik) + &
                                fact*overlap(1:nmatp,1:nmatp,ik) + &
                                (atot(1)/c)*pmat(1:nmatp,1:nmatp,1,ik) + &
                                (atot(2)/c)*pmat(1:nmatp,1:nmatp,2,ik) + &
                                (atot(3)/c)*pmat(1:nmatp,1:nmatp,3,ik)
    end do

#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
    
    call mt_h%release()

    if ( get_mathcalH ) call obtain_interstitial_contribution_mathcalH( &
      & first_kpt, last_kpt )

    if(tGen) then
      call timesec(timefinal)
      if(tDetail) tham = timefinal-timei
    end if

end subroutine UpdateHam
  !> Subroutine to calculate the interstitial contribution to `mathcalH` (used to
  !> obtain the forces on the ions in Ehrenfest Dynamics)
  subroutine obtain_interstitial_contribution_mathcalH( first_kpt, last_kpt )
    integer, intent(in) :: first_kpt, last_kpt

    integer     :: ik, is, ia, ias, ig, igl, j, ngp
    real(dp)    :: t1, t2, t3, t4, g(3)
    complex(dp) :: t5

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik,ngp,is,ia,ias,ig,igl,g,t1,t2,t3,t4,t5), &
!$OMP& SHARED(first_kpt,last_kpt,natoms), &
!$OMP& SHARED(nspecies,rmt,omega,atposc,idxas,apwalm,atot,pmat,ham_time), &
!$OMP& SHARED(ngk,nmat,vgkc,igkig,mathcalH), &
!$OMP& SHARED(nmatmax)
!$OMP DO
#endif    
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
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
  end subroutine obtain_interstitial_contribution_mathcalH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Subroutine to calculate the hamiltonian matrix for a given k-point
!> ik: the index of the k-point considered
!> nmatp:  the dimension of the matrix for this k-point (nmatp x nmatp)
!> calculate_mathcalH: if it is required to calculate the MT contributions 
!>                     to the auxiliary matrix mathcalH
  subroutine hamsetup( ik, nmatp, calculate_mathcalH )
    use constants, only: zzero, zone, zi
    use rttddft_GlobalVariables, only: mathcalH

    implicit none
    !> ik: the index of the k-point considered
    integer, intent(in)       :: ik
    !> `nmatp` is the dimension of the matrix for this `k-point` 
    !> (`nmatp` \( \times \) `nmatp`)
    integer, intent(in)       :: nmatp
    logical, intent(in)       :: calculate_mathcalH

    integer                   :: i, j, is, ia, ias, if3, ig, igl, io2
    integer                   :: j1, l3, m3, lm3, j3, maxnlo, maxaa
    integer                   :: iv(3)
    integer                   :: ngp
    real (dp)                 :: t1
    complex (dp)              :: zt
    complex (dp), allocatable :: hamcopy(:, :), aux(:, :)
    complex (dp), allocatable :: apwi(:, :), zm(:, :)

    if ( calculate_mathcalH ) mathcalH(:,:,:,:,ik) = zzero

    ! auxiliary variables
    ngp = ngk(1,ik)
    maxaa = mt_h%maxaa
    maxnlo = mt_h%maxnlo
    allocate( apwi(maxaa, ngp) )
    allocate( aux(nmatp, nmatp) )
    allocate( zm(maxaa, ngp) )
    allocate( hamcopy(nmatp, nmatp) )
    hamcopy(:,:) = zzero
    do is = 1, nspecies
      do ia = 1, natoms(is)
    ! APW-APW part
        ias = idxas (ia, is)
        apwi = zzero
        if3 = 0
        do l3 = 0, input%groundstate%lmaxmat
          do m3 = -l3, l3
          lm3 = idxlm (l3, m3)
            do io2 = 1, apword (l3, is)
              if3 = if3 + 1
              apwi(if3,1:ngp) = apwalm(1:ngp,io2,lm3,ias,ik)
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
              mathcalH(igl,ig,1:3,ias,ik) = &
                & zi*( vgkc(1:3,ig,1,ik)-vgkc(1:3,igl,1,ik) )*aux(igl,ig)
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
                mathcalH(j1:j3,ig,j,ias,ik) = zi*(vgkc(j,ig,1,ik))*aux(j1:j3,ig)
              end do
            end do
          end if
          do i = j1, j3
            hamcopy(1:ngp,i)=conjg(hamcopy(i,1:ngp))
            if ( calculate_mathcalH ) then
              do j = 1, 3
                mathcalH(1:ngp,i,j,ias,ik) = conjg(mathcalH(i,1:ngp,j,ias,ik))
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

    ham_time(1:nmatp,1:nmatp,ik) = hamcopy(1:nmatp,1:nmatp)
    
  end subroutine hamsetup

  !> Subroutine to calculate the overlap matrix for a given k-point.
  !> Based on `src/src_eigensystem/overlapsetup.f90`. 
  subroutine overlapsetup( ik, nmatp, calculate_mathcalB, calculate_mathcalH )
    use constants, only: zzero, zone, zi
    use physical_constants, only: alpha
    use rttddft_GlobalVariables, only: mathcalB, mathcalH, pmatmt

    implicit none

    !> ik: the index of the k-point considered
    integer, intent(in)       :: ik
    !> `nmatp` is the dimension of the matrix for this `k-point` 
    !> (`nmatp` \( \times \) `nmatp`)
    integer, intent(in)       :: nmatp
    !> calculate_mathcalB: if it is required to calculate the MT contributions 
    !>                     to the auxiliary matrix mathcalB
    !> calculate_mathcalH: if it is required to calculate the MT contributions 
    !>                     to the auxiliary matrix mathcalH
    logical, intent(in)       :: calculate_mathcalB
    logical, intent(in)       :: calculate_mathcalH

    integer                   :: i, is, ia, ias, if3, ig, j, j1, j2, igprime
    integer                   :: l, lm1, lm2, l3, m3, lm3
    integer                   :: io, io1, io2, maxaa, maxnlo, ilo, ilo1, ilo2
    integer                   :: ngp, iv(3)
    real (dp)                 :: t1
    real (dp), parameter      :: a2=0.5_dp*alpha**2
    complex (dp)              :: zt
    complex (dp), allocatable :: overlcopy(:, :), aux(:, :)
    complex (dp), allocatable :: apwi(:,:), zm(:,:), apwi2(:,:)

    if ( calculate_mathcalB ) mathcalB(:,:,:,:,ik) = -zi*pmatmt(:,:,:,:,ik)

    ngp = ngk(1,ik)
    maxaa = mt_h%maxaa
    maxnlo = mt_h%maxnlo
    allocate(apwi(maxaa,ngp))
    allocate(apwi2(ngp, maxaa) )
    allocate( overlcopy(nmatp, nmatp) )
    allocate( aux(nmatp, nmatp) )
    overlcopy(:,:) = zzero
    do is = 1, nspecies
      do ia = 1, natoms(is)
  ! APW-APW part
        ias = idxas (ia, is)
        apwi = zzero
        if3 = 0
        do l3 = 0, input%groundstate%lmaxmat
          do m3 = -l3, l3
          lm3 = idxlm (l3, m3)
            do io2 = 1, apword (l3, is)
              if3 = if3+1
              apwi(if3,1:ngp) = apwalm(1:ngp, io2, lm3, ias, ik)
            end do
          end do
        end do
        apwi2 = conjg( transpose( apwi ) )

        allocate( zm(ngp, maxaa) )
        zm(:,:) = zzero
        aux(:,:) = zzero
        if3 = 0
        do l3 = 0, input%groundstate%lmaxmat
          do m3 = -l3, l3
          lm3 = idxlm(l3, m3)
            do io2 = 1, apword(l3, is)
              do io1 = 1, apword (l3, is)
                zm(1:ngp,if3+io2) = zm(1:ngp,if3+io2) + &
                  & h1aa(io1,io2,l3,ias)*apwi2(1:ngp,if3+io1)
              enddo
              zm(1:ngp,if3+io2) = zm(1:ngp,if3+io2) + apwi2(1:ngp,if3+io2)
            end do
            if3 = if3 + apword(l3, is)
          end do
        end do
        call ZGEMM('C', 'C', ngp, ngp, maxaa, zone, apwi, maxaa, zm, ngp, &
          & zzero, aux, nmatp )
        overlcopy(1:ngp,1:ngp) = overlcopy(1:ngp,1:ngp) + aux(1:ngp,1:ngp)
        if ( calculate_mathcalB ) then
          do ig = 1, ngp
            do j = 1, 3
              mathcalB(1:ngp,ig,j,ias,ik) = mathcalB(1:ngp,ig,j,ias,ik) + &
                & zi*vgkc(j,ig,1,ik)*aux(1:ngp,ig)
            end do
          end do
        end if
        if ( calculate_mathcalH ) then
          do ig = 1, ngp
            do igprime = 1, ngp
              mathcalH(igprime,ig,1:3,ias,ik) = mathcalH(igprime,ig,1:3,ias,ik) + &
                & zi*( vgkc(1:3,ig,1,ik)-vgkc(1:3,igprime,1,ik) )*( &
                & fact*aux(igprime,ig) + (1/c)*( &
                & atot(1)*pmatmt(igprime,ig,1,ias,ik) + &
                & atot(2)*pmatmt(igprime,ig,2,ias,ik) + &
                & atot(3)*pmatmt(igprime,ig,3,ias,ik) ) )
            end do
          end do
        end if
        deallocate( zm )

  !What if it is, say, LAPW calculation without any local orbitals?
        if ( nlorb(is) /= 0 ) then
  ! APW-LO part
          !--Overlap--
          do ilo = 1, nlorb(is)
            l = lorbl(ilo, is)
            lm1 = idxlm(l,-l)
            lm2 = idxlm(l, l)
            j1 = ngp + idxlo(lm1, ilo, ias)
            j2 = ngp + idxlo(lm2, ilo, ias)
            aux(1:ngp,j1:j2) = zzero
            do io = 1, apword(l, is)
              aux(1:ngp,j1:j2) = aux(1:ngp,j1:j2) + &
                & conjg( apwalm(1:ngp,io,lm1:lm2,ias,ik) * &
                & ( oalo(io, ilo, ias) + h1loa(io, ilo, ias) ) )
            end do
            overlcopy(1:ngp,j1:j2) = overlcopy(1:ngp,j1:j2) + aux(1:ngp,j1:j2)
            if ( calculate_mathcalB ) then
              do ig = 1, ngp
                do j = 1, 3
                  mathcalB(j1:j2,ig,j,ias,ik) = mathcalB(j1:j2,ig,j,ias,ik) + &
                    & zi*vgkc(j,ig,1,ik)*conjg(aux(ig,j1:j2))
                end do
              end do
            end if
            if ( calculate_mathcalH ) then
              do ig = 1, ngp
                do j = j1, j2
                  mathcalH(ig,j,1:3,ias,ik) = mathcalH(ig,j,1:3,ias,ik) + &
                    & -zi*vgkc(1:3,ig,1,ik)*( aux(ig,j)*fact + (1/c)*( &
                    & atot(1)*pmatmt(ig,j,1,ias,ik) + &
                    & atot(2)*pmatmt(ig,j,2,ias,ik) + &
                    & atot(3)*pmatmt(ig,j,3,ias,ik) ) )
                  mathcalH(j,ig,1:3,ias,ik) = mathcalH(j,ig,1:3,ias,ik) + &
                    & zi*vgkc(1:3,ig,1,ik)*( conjg(aux(ig,j))*fact + (1/c)*( &
                    & atot(1)*pmatmt(j,ig,1,ias,ik) + &
                    & atot(2)*pmatmt(j,ig,2,ias,ik) + &
                    & atot(3)*pmatmt(j,ig,3,ias,ik) ) )
                end do
              end do
            end if
            do j = j1, j2
              overlcopy(j, 1:ngp)=conjg( overlcopy(1:ngp, j) )
            end do
          end do
  ! LO-LO part
          !--Overlap--
          do ilo1 = 1, nlorb(is)
            l = lorbl(ilo1,is)
            do ilo2 = 1, nlorb(is)
              if (lorbl(ilo2,is) .eq. l) Then
                lm1 = idxlm(l,-l)
                j1 = ngp + idxlo(lm1,ilo1,ias)
                j2 = ngp + idxlo(lm1,ilo2,ias)
                do lm2 = idxlm(l,-l), idxlm(l, l)
                    overlcopy(j1+lm2-lm1,j2+lm2-lm1) = &
                    & overlcopy(j1+lm2-lm1,j2+lm2-lm1) + &
                    & dcmplx(ololo(ilo1,ilo2,ias) + h1lolo(ilo1,ilo2,ias),0._dp)
                enddo
              end if
            end do ! do ilo2 = 1, nlorb(is)
          end do ! do ilo1 = 1, nlorb(is)
        endif ! if ( nlorb(is) /= 0 ) then
      end do ! do ia = 1, natoms(is)
    end do ! do is = 1, nspecies

  ! interstitial contributions
    if (input%groundstate%ValenceRelativity /= "none") then
      do j = 1, ngp
        do i = 1, j
          iv(:) = ivg(:,igkig(i,1,ik)) - ivg(:,igkig(j,1,ik))
          ig = ivgig(iv(1),iv(2),iv(3))
          if ((ig .gt. 0) .and. (ig .le. ngvec)) then
            t1 = 0.5_dp*dot_product(vgkc(:,i,1,ik),vgkc(:,j,1,ik))
            ! overlap
            zt = cfunig(ig)
            if ( input%groundstate%ValenceRelativity == 'iora*' ) then
              t1 = a2*t1
              zt = t1*m2effig(ig) + zt
            end if
            overlcopy(i,j) = overlcopy(i,j) + zt
            overlcopy(j,i) = conjg(overlcopy(i,j))
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
            overlcopy(i,j) = overlcopy(i,j) + cfunig(ig)
            overlcopy(j,i) = conjg(overlcopy(i,j))
          end if ! if ((ig .gt. 0) .and. (ig .le. ngvec))
        end do ! do i = 1, j
      end do ! do j = 1, ngp
    endif

    overlap(1:nmatp,1:nmatp,ik) = overlcopy(1:nmatp,1:nmatp)
    deallocate(apwi,apwi2)

  end subroutine overlapsetup


end module rttddft_HamiltonianOverlap
