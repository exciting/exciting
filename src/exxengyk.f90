
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
integer ngknr,ik,ist,jst
integer is,ia,ias,nrc,m,lmax
integer iv(3),iq,ig,igq0
real(8) cfq,v(3),t1
complex(8) zrho0,zt1
! automatic arrays
real(8) zn(nspecies)
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
! coefficient for long-range term
cfq=0.5d0*(omega/pi)**2
! set the nuclear charges to zero
zn(:)=0.d0
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(:,ikp),evalsvp)
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for occupied states for the input k-point
call genwfsv(.true.,ngk(1,ikp),igkig(:,1,ikp),evalsvp,apwalm,evecfv,evecsv, &
 wfmt1,wfir1)
! start loop over non-reduced k-point set
do ik=1,nkptnr
! generate G+k vectors
  call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
   tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-points
  call getevalsv(vklnr(:,ik),evalsvnr)
  call getevecfv(vklnr(:,ik),vgklnr,evecfv)
  call getevecsv(vklnr(:,ik),evecsv)
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
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vectors
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
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
  do jst=1,nstsv
    if (evalsvnr(jst).lt.efermi) then
      do ist=1,nstsv
        if (evalsvp(ist).lt.efermi) then
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(:,:,:,:,jst),wfmt1(:,:,:,:,ist), &
           wfir2(:,:,jst),wfir1(:,:,ist),zrhomt,zrhoir)
! calculate the Coulomb potential
          call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq, &
           sfacgq,zn,zrhomt,zrhoir,zvclmt,zvclir,zrho0) 
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
          t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
!$OMP CRITICAL
          evv=evv-0.5d0*occmax*wkpt(ikp)*(wkptnr(ik)*dble(zt1)+t1)
!$OMP END CRITICAL
! end loop over ist
        end if
      end do
! end loop over jst
    end if
  end do
! end loop over non-reduced k-point set
end do
!-----------------------------------------!
!    valence-core-valence contribution    !
!-----------------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do jst=1,spnst(is)
      if (spcore(jst,is)) then
        do m=-spk(jst,is),spk(jst,is)-1
! pass m-1/2 to wavefcr
          call wavefcr(lradstp,is,ia,jst,m,nrcmtmax,wfcr)
          do ist=1,nstsv
            if (evalsvp(ist).lt.efermi) then
! calculate the complex overlap density
              call vnlrhomt(.true.,is,wfcr(:,:,1),wfmt1(:,:,ias,1,ist), &
               zrhomt(:,:,ias))
              if (spinpol) then
                call vnlrhomt(.true.,is,wfcr(:,:,2),wfmt1(:,:,ias,2,ist),zfmt)
                zrhomt(:,1:nrc,ias)=zrhomt(:,1:nrc,ias)+zfmt(:,1:nrc)
              end if
! calculate the Coulomb potential
              call zpotclmt(ptnucl,lmaxvr,nrc,rcmt(:,is),0.d0,lmmaxvr, &
               zrhomt(:,:,ias),zvclmt(:,:,ias))
              zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
               zrhomt(:,:,ias),zvclmt(:,:,ias))
!$OMP CRITICAL
              ecv=ecv-occmax*wkpt(ikp)*dble(zt1)
!$OMP END CRITICAL
! end loop over ist
            end if
          end do
! end loop over m
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr)
deallocate(evalsvp,evalsvnr,evecfv,evecsv)
deallocate(sfacgknr,apwalm,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zfmt)
return
end subroutine

