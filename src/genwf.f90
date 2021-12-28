!
!
!
!BOP
! !ROUTINE: WFRelease
! !INTERFACE:
!
!
subroutine genWF(ik,wf)
! !USES:
! !DESCRIPTION:
! Generates wave functions in the WFType representation from first- and second-variational eigenvectors.
!
! !REVISION HISTORY:
!   Created 2021 (Andris)
!EOP
!BOC
!     use mod_kpoint, only : vkl
!     use mod_eigenvalue_occupancy, only : nstfv,nstsv
!     use mod_gkvector, only : ngk,vgkl,gkc,tpgkc,sfacgk,ngkmax
!     use mod_APW_LO, only : apwordmax
!     use mod_atoms, only : natmtot
!     use mod_muffin_tin, only : lmmaxapw
!     use modinput, only : input
 Use modinput
 Use mod_eigensystem
 use mod_kpoint
 use mod_eigenvalue_occupancy
 use mod_gkvector
 use mod_APW_LO
 use mod_atoms
 use mod_muffin_tin

 use modgw, only : kqset, Gkqset
 use mod_Gvector, only : ngrid, ngrtot, igfft
 Use mod_lattice, only : omega
 use constants, only : zzero, zone

!      use modmpi


implicit none
integer, intent (in) :: ik
type (WFType) :: wf

Complex (8), Allocatable :: apwalm (:, :, :, :)
Complex (8), Allocatable :: evecfv (:, :)
Complex (8), Allocatable :: evecsv (:, :)
integer :: ia,is,ias,if3,io,l,m,lm,wfsize,l1,l3,m1,m3,lm1,lm3,j1,j3
Complex (8), Allocatable :: apwi(:,:),wf1(:,:)
!     Complex (8), Allocatable :: zfft
real (8) :: t1
integer :: ifg,igk,j


Allocate(evecfv(nmatmax, nstfv))
Allocate(evecsv(nstsv, nstsv))
wfsize=wf%maxaa+wf%maxnlo
Allocate(wf1(wfsize,nstfv))
if (.not.allocated(wf%mt)) Allocate(wf%mt(wfsize,nstsv,natmtot))

if (.not.allocated(wf%gk)) Allocate(wf%gk(nmatmax, nstfv))
! get the eigenvalues/vectors from file for input k-point
!     Call getevalsv (vkl(:, ik), evalsv)
Call getevecfv (kqset%vkl(:, ik), Gkqset%vgkl(:, :, :, ik), evecfv)
!     Call getevecsv (vkl(:, ik), evecsv)
wf%gk = evecfv

! find the matching coefficients
Allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
Allocate(apwi(wf%maxaa,ngkmax))

!     Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), sfacgk(:, :, 1, ik), apwalm)
call match(Gkqset%ngk(1,ik), &
&          Gkqset%gkc(:,1,ik), &
&          Gkqset%tpgkc(:,:,1,ik), &
&          Gkqset%sfacgk(:,:,1,ik),&
&          apwalm)


Do is = 1, nspecies
!       n = lmmaxvr * nrcmt (is)
  Do ia = 1, natoms (is)
    ias = idxas (ia, is)
    wf1=zzero

    if3=0
    Do l = 0, input%groundstate%lmaxmat
      Do io = 1, apword (l, is)
        Do m = - l, l
          lm = idxlm (l, m)
          if3=if3+1
          apwi(if3,1:Gkqset%ngk(1, ik))=apwalm(1:Gkqset%ngk(1, ik), io, lm, ias)
        End Do
      End Do
    End Do

! APW coefficients
     call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                'N', &           ! TRANSB = 'N'  op( B ) = B.
                 if3, &          ! M ... rows of op( A ) = rows of C
                 nstfv, &           ! N ... cols of op( B ) = cols of C
                 Gkqset%ngk(1,ik), &        ! K ... cols of op( A ) = rows of op( B )
                 zone, &          ! alpha
                 apwi, &        ! A
                 wf%maxaa,&           ! LDA ... leading dimension of A
                 evecfv, &           ! B
                 nmatmax, &          ! LDB ... leading dimension of B
                 zzero, &          ! beta
                 wf1, &  ! C
                 wfsize &      ! LDC ... leading dimension of C
                 )
! LO coefficients
  if (wf%losize(is).gt.0) then
    l1 = lorbl (1, is)
    l3 = lorbl (nlorb(is), is)
    lm1=idxlm (l1,-l1)
    lm3=idxlm (l3,l3)
    j1= Gkqset%ngk(1, ik) + idxlo (lm1, 1, ias)
    j3= Gkqset%ngk(1, ik) + idxlo (lm3, nlorb(is), ias)
    wf1(if3+1:if3+wf%losize(is),1:nstfv)=evecfv(j1:j3,1:nstfv)
  endif
  ! wf%mt(wfsize,nstsv,natmtot)
  wf%mt(:,:,ias)=wf1(:,:)

! Apply the second variation
if (.false.) then
     call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                'N', &           ! TRANSB = 'N'  op( B ) = B.
                 wfsize, &          ! M ... rows of op( A ) = rows of C
                 nstsv, &           ! N ... cols of op( B ) = cols of C
                 nstfv, &        ! K ... cols of op( A ) = rows of op( B )
                 zone, &          ! alpha
                 wf1, &        ! A
                 wfsize,&           ! LDA ... leading dimension of A
                 evecsv(1,1), &           ! B
                 nstsv, &          ! LDB ... leading dimension of B
                 zzero, &          ! beta
                 wf%mt(1,1,ias), &  ! C
                 wfsize &      ! LDC ... leading dimension of C
                 )
endif

  End Do
End Do

Deallocate(apwi)
Deallocate(apwalm)
Deallocate(wf1)

!if (.not.allocated(wf%ir)) Allocate(wf%ir(ngrtot,nstsv))
!
!!     interstitial part
! t1 = 1 / sqrt(omega)
! Do j = 1, nstsv
!   wf%ir(:,j) = 0.d0
!! spin-unpolarised wavefunction
!   Do igk = 1, Gkqset%ngk (1, ik)
!     ifg = igfft (Gkqset%igkig(igk, 1, ik))
!     wf%ir(ifg, j) = t1*evecfv(igk, j)
!   End Do
!! Fourier transform wavefunction to real-space
!   Call zfftifc (3, ngrid, 1, wf%ir(:, j))
! End Do
!! wf%ir(:,:) = 0.d0
!! wf%mt(:,:,:)=0.d0
Deallocate(evecfv)
Deallocate(evecsv)

end subroutine genWF

