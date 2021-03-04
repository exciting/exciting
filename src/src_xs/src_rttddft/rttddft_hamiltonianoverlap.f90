! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created Apr 2019 (Ronaldo)
! Reference: https://arxiv.org/abs/2102.02630

!> Module that manages the hamiltonian and overlap in RT-TDDFT
!> TODO(Ronaldo): Refactor to reduce the number of global variables
module rtddft_HamiltonianOverlap

  use modmpi
  use modinput, only: input
  use mod_APW_LO, only: apword, nlorb, lorbl
  use mod_kpoint, only: nkpt
  use mod_gkvector, only: ngk, vgkc, igkig
  use mod_potential_and_density, only: veffig, meffig, m2effig
  use mod_atoms, only: nspecies, natoms, idxas, natmtot, atposc
  use mod_gvector, only: ivg, ivgig, cfunig, ngvec
  use mod_muffin_tin, only: idxlm, rmt
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_lattice, only: omega
  use mod_eigensystem, only: nmat, nmatmax, hloloij, idxlo, h1aa, h1loa, h1lolo, oalo, ololo
  use rttddft_GlobalVariables, only: ham_time, ham_past, overlap, apwalm, &
    & atot, pmat, evecfv_time, timesecRTTDDFT
  use precision, only: dp
  use physical_constants, only: c
  use mod_eigensystem, only: MTHamiltonianList, MTInitAll, MTNullify

  type (MTHamiltonianList) :: mt_h

  private

  public :: UpdateHam

contains

  !> In UpdateHam, we obtain the hamiltonian (and if requested, the overlap) at time \( t \).
  !> The subroutines here are based on src\src_eigensystem\hamiltonsetup.f90
  !> and src\src_eigensystem\overlapsetup.f90
  !> We decided not to use them because we would need to initialize
  !> the variable "Type (evsystem) :: system" at all steps, which, repeated so
  !> many times over the propagation of the wavefunctions, would include
  !> undesirable delays.
  !> @param[in]   predcorr            tells if we are in the loop of the
  !>                                  predictor-Corrector scheme
  !> @param[in]   calculateOverlap    tells if we need to calculate the overlap
  !> @param[in]   timeGen             tells if we want a general timing
  !> @param[in]   timeDetail          tells if we want a detailed timing
  !>                                  (only works if timeGen is true)
  !> @param[in]   timeini             time (in seconds) elapsed since exciting
  !>                                  was started (employed for timing)
  !> @param[out]   timefinal          time (in seconds) after executing this subroutine
  !> @param[out]   thmlint            time spent to execute hmlint
  !> @param[out]   tham               time spent after executing hmlint until
  !>                                  the end of this subroutine
  subroutine UpdateHam( predcorr, calculateOverlap, &
    & timeGen, timeDetail, timeini, timefinal, thmlint, tham )

    implicit none

    !> predcorr: tells if we are in the loop of the predictor-Corrector scheme
    !> calculateOverlap:  tells if we need to calculate the overlap
    logical, intent(in)               :: predcorr, calculateOverlap
    !> timeGen: tells if we want a general timing
    !> timeDetail: tells if we want a detailed timing (only works if timeGen is true)
    logical, intent(in), optional     :: timeGen, timeDetail
    !> timeini: time (in seconds) elapsed since exciting was started (employed for timing)
    real(dp), intent(in), optional    :: timeini
    !> timefinal: time (in seconds) after executing this subroutine
    !> thmlint  : time spent to execute hmlint
    !> tham     : time spent after executing hmlint until the end of this subroutine
    real(dp), intent(out), optional   :: timefinal, thmlint, tham

    integer               :: ik, nmatp, first_kpt, last_kpt
    real(dp)              :: timei, timef, fact
    logical               :: tGen,tDetail

  ! factor that multiplies the overlap matrix (when we compute the hamiltonian)
  fact = dot_product(atot, atot) / (2._dp * c**2)

  ! Check optional arguments
  if ( present(timeGen) ) then
    tGen = timeGen
    if ( present(timeDetail) ) then
      tDetail = timeDetail
    else
      tDetail = .false.
    end if
  else
    tGen = .False.
    tDetail = .False.
  end if
  if( tGen ) timei = timeini

  !if( tDetail ) call timesecRTTDDFT(timei,timef,tgenpmatbasis)
  call MTNullify(mt_h)
  call MTInitAll(mt_h)
  call hmlint(mt_h)

  if ( tDetail ) call timesecRTTDDFT( timei, timef, thmlint )

  if ( .not. predcorr ) ham_past(:,:,:) = ham_time(:,:,:)

#ifdef MPI
  first_kpt = firstk(rank)
  last_kpt = lastk(rank)
#else
  first_kpt = 1
  last_kpt = nkpt
#endif

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik,nmatp), &
!$OMP& SHARED(first_kpt,last_kpt,calculateOverlap,nkpt,natmtot,natoms), &
!$OMP& SHARED(nspecies,rmt,omega,atposc,idxas,apwalm,fact,atot,pmat,ham_time), &
!$OMP& SHARED(overlap,ngk,nmat,vgkc,igkig,input), &
!$OMP& SHARED(nmatmax)
!$OMP DO
#endif
  do ik = first_kpt, last_kpt
    nmatp = nmat(1,ik)
    call hamsetup(ik,nmatp)
    if ( calculateOverlap ) call overlapsetup(ik,nmatp)

    ! Include the part of the vector potential in the hamiltonian
    ham_time(1:nmatp,1:nmatp,ik) = ham_time(1:nmatp,1:nmatp,ik) + &
                                fact*overlap(1:nmatp,1:nmatp,ik) + &
                                (atot(1)/c)*pmat(1:nmatp,1:nmatp,1,ik) + &
                                (atot(2)/c)*pmat(1:nmatp,1:nmatp,2,ik) + &
                                (atot(3)/c)*pmat(1:nmatp,1:nmatp,3,ik)
  end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  if(tGen) then
    call timesec(timefinal)
    if(tDetail) tham = timefinal-timei
  end if

end subroutine updateham
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Subroutine to calculate the hamiltonian matrix for a given k-point
!> @param[in]   ik      the index of the k-point considered
!> @param[in]   nmatp   the dimension of the matrix for this k-point (nmatp x nmatp)
  subroutine hamsetup(ik,nmatp)
    use constants, only: zzero, zone
    implicit none
    !> ik: the index of the k-point considered
    !> nmatp:  the dimension of the matrix for this k-point (nmatp x nmatp)
    integer, intent(in)       :: ik, nmatp

    integer                   :: i, j, is, ia, ias, if3, ig, io2
    integer                   :: j1, l3, m3, lm3, j3, maxnlo, maxaa
    integer                   :: iv(3)
    integer                   :: ngp
    real (dp)                 :: t1
    complex (dp)              :: zt
    complex (dp), allocatable :: hamcopy(:, :)
    complex (dp), allocatable :: apwi(:, :), zm(:, :)

    ! auxiliary variables
    ngp = ngk(1,ik)
    maxaa = mt_h%maxaa
    maxnlo = mt_h%maxnlo
    allocate( apwi(maxaa, ngp) )
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
          & zone, hamcopy, nmatp )

    !What if it is, say, LAPW calculation without any local orbitals?
        if ( nlorb(is) /= 0 ) then
    ! APW-LO part
          l3 = lorbl(1,is)
          lm3 = idxlm(l3,-l3)
          call ZGEMM( 'N', 'N', mt_h%losize(is), ngp, maxaa, zone, &
            & mt_h%main%loa(:,:,ias), maxnlo, apwi, maxaa, zone, &
            & hamcopy(ngp+idxlo(lm3,1,ias),1), nmatp )
          j1 = ngp + idxlo( lm3, 1, ias )
          j3 = j1 + mt_h%losize(is) - 1
          do i = j1, j3
            hamcopy(1:ngp,i)=conjg(hamcopy(i,1:ngp))
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine to calculate the overlap matrix for a given k-point
  !> @param[in]   ik      the index of the k-point considered
  !> @param[in]   nmatp   the dimension of the matrix for this k-point (nmatp x nmatp)
  subroutine overlapsetup( ik, nmatp )
    use constants, only: zzero, zone
    use physical_constants, only: alpha

    implicit none
    !> ik: the index of the k-point considered
    !> nmatp:  the dimension of the matrix for this k-point (nmatp x nmatp)
    integer, intent(in)       :: ik, nmatp

    integer                   :: i, is, ia, ias, if3, ig, j, j1, j2
    integer                   :: l, lm1, lm2, l3, m3, lm3
    integer                   :: io, io1, io2, maxaa, maxnlo, ilo, ilo1, ilo2
    integer                   :: ngp, iv(3)
    real (dp)                 :: t1
    real (dp), parameter      :: a2=0.5_dp*alpha**2
    complex (dp)              :: zt
    complex (dp), allocatable :: overlcopy(:, :)
    complex (dp), allocatable :: apwi(:,:), zm(:,:), apwi2(:,:)

    ngp = ngk(1,ik)
    maxaa = mt_h%maxaa
    maxnlo = mt_h%maxnlo
    allocate(apwi(maxaa,ngp))
    allocate(apwi2(ngp,maxaa))
    allocate(overlcopy(nmatp,nmatp))
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
          & zone, overlcopy, nmatp )
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
            do io = 1, apword(l, is)
              overlcopy(1:ngp,j1:j2) = overlcopy(1:ngp,j1:j2) + &
                & conjg( apwalm(1:ngp,io,lm1:lm2,ias,ik) * &
                & ( oalo(io, ilo, ias) + h1loa(io, ilo, ias) ) )
            end do
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


end module rtddft_HamiltonianOverlap
