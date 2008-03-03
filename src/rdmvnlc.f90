
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmvnlc(ikp,vnl)
! calculate non-local matrix elements for minimisation w.r.t. evecsv
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vnl(nstsv,nstsv,nstsv,nkptnr)
! local variables
integer ngknr,ik,ist1,ist2,ist3
integer lmax,ig,iq,igq0,iv(3)
real(8) cfq,v(3),t1
complex(8) zrho0,zt1,zt2,zrho1
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
real(8), allocatable :: evalsvl(:)
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
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zpchg(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
! external functions
complex(8) zfinp
external zfinp
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
allocate(evalsvl(nstsv))
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
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zpchg(natmtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! set the point charges to zero
zpchg(:)=0.d0
! get the eigenvectors and values from file
call getevalsv(vkl(1,ikp),evalsvl)
call getevecfv(vkl(1,ikp),vgkl(1,1,ikp,1),evecfv)
call getevecsv(vkl(1,ikp),evecsv)
! find the matching coefficients
call match(ngk(ikp,1),gkc(1,ikp,1),tpgkc(1,1,ikp,1),sfacgk(1,1,ikp,1),apwalm)
! calculate the wavefunctions for all states for the input k-point
call genwfsv(.false.,ngk(ikp,1),igkig(1,ikp,1),evalsvl,apwalm,evecfv,evecsv, &
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
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
! compute the required spherical Bessel functions
  lmax=lmaxvr+npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
! calculate the wavefunctions for all states for non-reduced k-point ik
  call genwfsv(.false.,ngknr,igkignr,evalsvnr,apwalm,evecfv,evecsv,wfmt2,wfir2)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
  do ist1=1,nstsv
    do ist2=1,nstsv
! calculate the complex overlap density
      call vnlrho(.true.,wfmt2(1,1,1,1,ist2),wfmt1(1,1,1,1,ist1), &
       wfir2(1,1,ist2),wfir1(1,1,ist1),zrhomt,zrhoir)
! compute the potential and G=0 coefficient of the density
      call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
       zpchg,zrhomt,zrhoir,zvclmt,zvclir,zrho0)
      zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
      t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
      vnl(ist1,ist1,ist2,ik)=wkptnr(ik)*dble(zt1)+t1
      do ist3=1,nstsv
        if (ist1.gt.ist3) then
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(1,1,1,1,ist2),wfmt1(1,1,1,1,ist3), &
           wfir2(1,1,ist2),wfir1(1,1,ist3),zrhomt,zrhoir)
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
! compute the density coefficient of the smallest G+q-vector
          call zrhoqint(gqc(igq0),ylmgq(1,igq0),ngvec,sfacgq(igq0,1),zrhomt, &
           zrhoir,zrho1)
          zt2=cfq*wiq2(iq)*(conjg(zrho1)*zrho0)
          vnl(ist3,ist1,ist2,ik)=wkptnr(ik)*zt1+zt2
! end loop over ist3
        end if
      end do
! end loop over ist2
    end do
! end loop over ist1
  end do
! calculate the lower diagonal
  do ist1=1,nstsv
    do ist3=1,nstsv
      if (ist1.lt.ist3) then
        vnl(ist3,ist1,:,ik)=conjg(vnl(ist1,ist3,:,ik))
      end if
    end do
  end do
! end loop over non-reduced k-point set
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr,vgqc,tpgqc,gqc,jlgqr)
deallocate(evalsvl,evalsvnr)
deallocate(apwalm,evecfv,evecsv,sfacgknr,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zpchg,zvclmt,zvclir)
return
end subroutine
!EOC

