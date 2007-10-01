
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengyk(ikp,evv,ecv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(inout) :: evv
real(8), intent(inout) :: ecv
! local variables
integer ngknr,ik,ist1,ist2
integer is,ia,ias,nr,m,lmax
integer iv(3),iq,ig,igq0
real(8) occ,cfq,v(3),t1
complex(8) zrho0,zt1
! allocatable arrays
integer, allocatable :: igkignr(:)
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
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: wfcr(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zpchg(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
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
allocate(jlgqr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(evalsvp(nstsv))
allocate(evalsvnr(nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(ylmgq(lmmaxvr,ngvec))
allocate(sfacgq(ngvec,natmtot))
allocate(zpchg(natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
allocate(wfcr(lmmaxvr,nrcmtmax,2))
allocate(zfmt(lmmaxvr,nrcmtmax))
if (spinpol) then
  occ=1.d0
else
  occ=2.d0
end if
! coefficient for long-range term
cfq=0.5d0*(omega/pi)**2
! set the point charges to zero
zpchg(:)=0.d0
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(1,ikp),evalsvp)
call getevecfv(vkl(1,ikp),vgkl(1,1,ikp,1),evecfv)
call getevecsv(vkl(1,ikp),evecsv)
! find the matching coefficients
call match(ngk(ikp,1),gkc(1,ikp,1),tpgkc(1,1,ikp,1),sfacgk(1,1,ikp,1),apwalm)
! calculate the wavefunctions for occupied states for the input k-point
call genwfsv(.true.,ngk(ikp,1),igkig(1,ikp,1),evalsvp,apwalm,evecfv,evecsv, &
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
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkcnr(:,ik)
  do ig=1,ngvec
! determine G+q vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(1,ig),gqc(ig),tpgqc(1,ig))
! spherical harmonics for G+q-vectors
    call genylm(lmaxvr,tpgqc(1,ig),ylmgq(1,ig))
  end do
! structure factor for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
! compute the required spherical Bessel functions
  lmax=lmaxvr+npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
! calculate the wavefunctions for occupied states
  call genwfsv(.true.,ngknr,igkignr,evalsvnr,apwalm,evecfv,evecsv,wfmt2,wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
  do ist2=1,nstsv
    if (evalsvnr(ist2).lt.efermi) then
      do ist1=1,nstsv
        if (evalsvp(ist1).lt.efermi) then
! calculate the complex overlap density
          call vnlrho(wfmt2(1,1,1,1,ist2),wfmt1(1,1,1,1,ist1),wfir2(1,1,ist2), &
           wfir1(1,1,ist1),zrhomt,zrhoir)
! calculate the Coulomb potential
          call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq, &
           sfacgq,zpchg,zrhomt,zrhoir,zvclmt,zvclir,zrho0) 
          zt1=zfinp(zrhomt,zvclmt,zrhoir,zvclir)
          t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
!$OMP CRITICAL
          evv=evv-0.5d0*occ*wkpt(ikp)*(wkptnr(ik)*dble(zt1)+t1)
!$OMP END CRITICAL
! end loop over ist1
        end if
      end do
! end loop over ist2
    end if
  end do
! end loop over non-reduced k-point set
end do
!-----------------------------------------!
!    valence-core-valence contribution    !
!-----------------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nr=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist2=1,spnst(is)
      if (spcore(ist2,is)) then
        do m=-spk(ist2,is),spk(ist2,is)-1
! pass m-1/2 to wavefcr
          call wavefcr(lradstp,is,ia,ist2,m,nrcmtmax,wfcr)
          do ist1=1,nstsv
            if (evalsvp(ist1).lt.efermi) then
! calculate the complex overlap density
              call vnlrhomt(is,wfcr(1,1,1),wfmt1(1,1,ias,1,ist1), &
               zrhomt(1,1,ias))
              if (spinpol) then
                call vnlrhomt(is,wfcr(1,1,2),wfmt1(1,1,ias,2,ist1),zfmt)
                zrhomt(:,1:nr,ias)=zrhomt(:,1:nr,ias)+zfmt(:,1:nr)
              end if
! calculate the Coulomb potential
              call zpotclmt(lmaxvr,nr,rcmt(1,is),zpchg(ias),lmmaxvr, &
               zrhomt(1,1,ias),zvclmt(1,1,ias))
              zt1=zfmtinp(lmaxvr,nr,rcmt(1,is),lmmaxvr,zrhomt(1,1,ias), &
               zvclmt(1,1,ias))
!$OMP CRITICAL
              ecv=ecv-occ*wkpt(ikp)*dble(zt1)
!$OMP END CRITICAL
! end loop over ist1
            end if
          end do
! end loop over m
        end do
! end loop over ist2
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr,vgqc,tpgqc,gqc,jlgqr)
deallocate(evalsvp,evalsvnr,evecfv,evecsv)
deallocate(sfacgknr,apwalm,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr)
deallocate(zrhomt,zrhoir,zpchg,zvclmt,zvclir,zfmt)
return
end subroutine

