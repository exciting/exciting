
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnhf(ikp,evecsvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(inout) :: evecsvp(nstsv,nstsv)
! local variables
integer is,ia,ias,ir,irc
integer ngknr,ik,jk,ist1,ist2,ist3
integer iq,ig,iv(3),igq0
integer lmax,lwork,info
real(8) cfq,v(3),t1
complex(8) zrho01,zrho02,zt1,zt2
! automatic arrays
real(8) zn(nspecies)
complex(8) sfacgq0(natmtot)
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
real(8), allocatable :: jlgq0r(:,:,:)
real(8), allocatable :: evalsvp(:)
real(8), allocatable :: evalsvnr(:)
real(8), allocatable :: rwork(:)
real(8), allocatable :: rfmt(:,:,:)
complex(8), allocatable :: h(:,:)
complex(8), allocatable :: vmat(:,:)
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
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
complex(8), allocatable :: zfmt(:,:)
complex(8), allocatable :: work(:)
! external functions
complex(8) zfinp
external zfinp
!$OMP CRITICAL
write(*,'("Info(seceqnhf): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL
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
allocate(jlgq0r(0:lmaxvr,nrcmtmax,nspecies))
allocate(evalsvp(nstsv))
allocate(evalsvnr(nstsv))
allocate(rwork(3*nstsv))
allocate(h(nstsv,nstsv))
allocate(vmat(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(ylmgq(lmmaxvr,ngvec))
allocate(sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
allocate(zfmt(lmmaxvr,nrcmtmax))
lwork=2*nstsv
allocate(work(lwork))
allocate(rfmt(lmmaxvr,nrcmtmax,natmtot))
! coefficient of long-range term
cfq=0.5d0*(omega/pi)**2
! set the point charges to zero
zn(:)=0.d0
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(:,ikp),evalsvp)
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for all states for the input k-point
call genwfsv(.false.,ngk(1,ikp),igkig(:,1,ikp),evalsvp,apwalm,evecfv,evecsvp, &
 wfmt1,wfir1)
! compute the new kinetic matrix elements
call zgemm('N','N',nstsv,nstsv,nstsv,zone,kinmatc(:,:,ikp),nstsv,evecsvp, &
 nstsv,zzero,vmat,nstsv)
call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvp,nstsv,vmat,nstsv,zzero,h, &
 nstsv)
! convert muffin-tin Coulomb potential to spherical coordinates
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,vclmt(:,ir,ias),1, &
       0.d0,rfmt(:,irc,ias),1)
    end do
  end do
end do
! compute the Coulomb matrix elements and add
call genvmatk(rfmt,vclir,wfmt1,wfir1,vmat)
h(:,:)=h(:,:)+vmat(:,:)
! zero the non-local matrix elements for passed k-point
vmat(:,:)=0.d0
! start loop over non-reduced k-point set
do ik=1,nkptnr
! find the equivalent reduced k-point
  iv(:)=ivknr(:,ik)
  jk=ikmap(iv(1),iv(2),iv(3))
! generate the G+k vectors
  call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
   tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-point
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
! spherical harmonics for G+q-vector
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  lmax=lmaxvr+npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
  call genjlgq0r(gqc(igq0),jlgq0r)
! calculate the wavefunctions for all states
  call genwfsv(.false.,ngknr,igkignr,evalsvnr,apwalm,evecfv,evecsv,wfmt2,wfir2)
  do ist3=1,nstsv
    if (occsv(ist3,jk).gt.epsocc) then
      do ist2=1,nstsv
! calculate the complex overlap density
        call vnlrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist2), &
         wfir2(:,:,ist3),wfir1(:,:,ist2),zrhomt,zrhoir)
! calculate the Coulomb potential
        call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq, &
         sfacgq,zn,zrhomt,zrhoir,zvclmt,zvclir,zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
        do ist1=1,ist2
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist1), &
           wfir2(:,:,ist3),wfir1(:,:,ist1),zrhomt,zrhoir)
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
! compute the density coefficient of the smallest G+q-vector
          call zrhogp(gqc(igq0),jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt,zrhoir, &
           zrho01)
          zt2=cfq*wiq2(iq)*(conjg(zrho01)*zrho02)
          t1=occsv(ist3,jk)/occmax
          vmat(ist1,ist2)=vmat(ist1,ist2)-t1*(wkptnr(ik)*zt1+zt2)
        end do
      end do
    end if
  end do
! end loop over non-reduced k-point set
end do
! add the non-local matrix elements to Hamiltonian
h(:,:)=h(:,:)+vmat(:,:)
! diagonalise the Hartree-Fock Hamiltonian (eigenvalues in global array)
call zheev('V','U',nstsv,h,nstsv,evalsv(:,ikp),work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnhf): diagonalisation of the Hartree-Fock Hamiltonian &
   &failed")')
  write(*,'(" for k-point ",I8)') ikp
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! apply unitary transformation to second-variational states
evecsv(:,:)=evecsvp(:,:)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,h,nstsv,zzero,evecsvp, &
 nstsv)
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(evalsvp,evalsvnr,evecfv,evecsv,rwork)
deallocate(h,vmat,apwalm,sfacgknr,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zfmt,work,rfmt)
return
end subroutine

