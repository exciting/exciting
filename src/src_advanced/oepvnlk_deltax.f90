
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!Created March 2008 (DIN)
!Revision March 2014 (UW)


subroutine oepvnlk_deltax(ikp,vnlvv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vnlvv(nstsv)
! local variables
integer ngknr,ik,ist2,ist3
integer is,ia,ias,ic,m1,m2,lmax
integer nr,iq,ig,iv(3),igq0
real(8) v(3),cfq
complex(8) zrho01,zrho02,zt1,zt2
Complex (8) sfacgq0 (natmtot)
! allocatable arrays
integer, allocatable :: igkignr(:)
 Real (8), Allocatable :: jlgq0r (:, :, :)
real(8), allocatable :: vgklnr(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:)
real(8), allocatable :: tpgqc(:,:)
real(8), allocatable :: gqc(:)
real(8), allocatable :: jlgqr(:,:,:)
real(8), allocatable :: evalsvp(:)
real(8), allocatable :: evalsvnr(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: ylmgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: wfcr1(:,:,:)
complex(8), allocatable :: wfcr2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zpchg(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
complex(8), allocatable :: zvcltp(:,:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax))
allocate(vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax))
allocate(tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec))
allocate(tpgqc(2,ngvec))
allocate(gqc(ngvec))
allocate(jlgqr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1,ngvec,nspecies))
Allocate (jlgq0r(0:input%groundstate%lmaxvr, nrcmtmax, nspecies))
allocate(evalsvp(nstsv))
allocate(evalsvnr(nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec))
allocate(sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(wfcr1(lmmaxvr,nrcmtmax,2))
allocate(wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zpchg(natmtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
allocate(zvcltp(lmmaxvr,nrcmtmax))
allocate(zfmt(lmmaxvr,nrcmtmax))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! set the point charges to zero
zpchg(:)=0.d0
vnlvv(:)=0.d0
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(1,ikp),evalsvp)
call getevecfv(vkl(1,ikp),vgkl(1,1,ikp,1),evecfv)
call getevecsv(vkl(1,ikp),evecsv)
! find the matching coefficients
call match(ngk(ikp,1),gkc(1,ikp,1),tpgkc(1,1,ikp,1),sfacgk(1,1,ikp,1),apwalm)
! calculate the wavefunctions for all states for the input k-point
call genwfsv(.false.,ngk(ikp,1),igkig(1,ikp,1),evalsvp,apwalm,evecfv,evecsv, &
 wfmt1,wfir1)
! start loop over non-reduced k-point set
do ik=1,nkptnr
! generate G+k vectors
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
   tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-points
  call getevalsv(vklnr(1,ik),evalsvnr)
  call getevecfv(vklnr(1,ik),vgklnr,evecfv)
  call getevecsv(vklnr(1,ik),evecsv)
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! determine q-vector
  iv(:)=ivk(:,ikp)-ivknr(:,ik)
  iv(:)=modulo(iv(:),input%groundstate%ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkcnr(:,ik)
  do ig=1,ngvec
! determine G+q vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(1,ig),gqc(ig),tpgqc(1,ig))
! spherical harmonics for G+q-vector
    call genylm(input%groundstate%lmaxvr,tpgqc(1,ig),ylmgq(1,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0 (:) = sfacgq (igq0, :)

! compute the required spherical Bessel functions
  lmax=input%groundstate%lmaxvr+input%groundstate%npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
  Call genjlgq0r (gqc(igq0), jlgq0r)
! calculate the wavefunctions for occupied states
  call genwfsv(.true.,ngknr,igkignr,evalsvnr,apwalm,evecfv,evecsv,wfmt2,wfir2)
  do ist3=1,nstsv
    if (evalsvnr(ist3).lt.efermi) then
      do ist2=1,nstsv
!        if (evalsvp(ist2).gt.efermi) then
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(1,1,1,1,ist3),wfmt1(1,1,1,1,ist2), &
           wfir2(1,1,ist3),wfir1(1,1,ist2),zrhomt,zrhoir)
! calculate the Coulomb potential
          call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq, &
           sfacgq,zpchg,zrhomt,zrhoir,zvclmt,zvclir,zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
! ## try: replace following lines by diagonal elements ##
!          do ist1=1,nstsv
!            if (evalsvp(ist1).lt.efermi) then
! calculate the complex overlap density
              call vnlrho(.true.,wfmt2(1,1,1,1,ist3),wfmt1(1,1,1,1,ist2), &
               wfir2(1,1,ist3),wfir1(1,1,ist2),zrhomt,zrhoir)
              zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
! compute the density coefficient of the smallest G+q-vector
                Call zrhogp (gqc(igq0), jlgq0r, ylmgq(:,igq0), sfacgq0, &
                  &  zrhomt, zrhoir, zrho01)

              zt2=cfq*wiq2(iq)*(conjg(zrho01)*zrho02)
              vnlvv(ist2)=vnlvv(ist2)-(wkptnr(ik)*zt1+zt2)
!            end if
!          end do

! end loop over ist2
!       end if
      end do
! end loop over ist3
    end if
  end do
! end loop over non-reduced k-point set
end do
! begin loops over atoms and species
do is=1,nspecies
  nr=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,spnst(is)
      if (spcore(ist3,is)) then
        do m1=-spk(ist3,is),spk(ist3,is)-1
! pass m-1/2 to wavefcr
          call wavefcr(input%groundstate%lradstep,is,ia,ist3,m1,nrcmtmax,wfcr1)
          do ist2=1,nstsv
!            if (evalsvp(ist2).gt.efermi) then
! calculate the complex overlap density
              call vnlrhomt(.true.,is,wfcr1(1,1,1),wfmt1(1,1,ias,1,ist2), &
               zrhomt(1,1,ias))
              !if (spinpol) then
               If (associated(input%groundstate%spin)) Then
                call vnlrhomt(.true.,is,wfcr1(1,1,2),wfmt1(1,1,ias,2,ist2),zfmt)
                zrhomt(:,1:nr,ias)=zrhomt(:,1:nr,ias)+zfmt(:,1:nr)
              end if
! calculate the Coulomb potential
              Call zpotclmt (input%groundstate%ptnucl, &
             &  input%groundstate%lmaxvr, nr, rcmt(:, is), &
             &  0.d0, lmmaxvr, zrhomt(:, :, ias), zfmt)
! convert muffin-tin potential to spherical coordinates
              call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtapw, &
               lmmaxapw,zfmt,lmmaxvr,zzero,zvcltp,lmmaxvr)
!-------------------------------------------!
!     valence-core-valence contribution     !
!-------------------------------------------!
!              do ist1=1,nstsv
!                if (evalsvp(ist1).lt.efermi) then
! calculate the complex overlap density
                  call vnlrhomt(.false.,is,wfcr1(1,1,1),wfmt1(1,1,ias,1,ist2), &
                   zrhomt(1,1,ias))
                  If (associated(input%groundstate%spin)) Then
                  !if (spinpol) then
                    call vnlrhomt(.false.,is,wfcr1(1,1,2), &
                     wfmt1(1,1,ias,2,ist2),zfmt)
                    zrhomt(:,1:nr,ias)=zrhomt(:,1:nr,ias)+zfmt(:,1:nr)
                  end if
                  zt1=zfmtinp(.false.,input%groundstate%lmaxvr,nr,rcmt(1,is),lmmaxvr, &
                   zrhomt(1,1,ias),zvcltp)
                  vnlvv(ist2)=vnlvv(ist2)-zt1
!                end if
!              end do
! end loop over ist2
!            end if
          end do
! end loops over ist3 and m1
        end do
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr,vgqc,tpgqc,gqc,jlgqr)
deallocate(evalsvp,evalsvnr,evecfv,evecsv)
deallocate(apwalm,sfacgknr,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr1,wfcr2)
deallocate(zrhomt,zrhoir,zpchg,zvclmt,zvclir,zvcltp,zfmt)
return
end subroutine
!EOC

