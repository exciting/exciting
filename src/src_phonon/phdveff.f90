
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdveff(iph,iq,veffmtp,veffirp,dveffmt,dveffir)
use modmain
implicit none
! arguments
integer, intent(in) :: iph
integer, intent(in) :: iq
real(8), intent(in) :: veffmtp(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: veffirp(ngrtot)
complex(8), intent(out) :: dveffmt(lmmaxvr,nrcmtmax,natmtot0)
complex(8), intent(out) :: dveffir(ngrtot0)
! local variables
integer is,ia,ja,ias,jas
integer ir,irc,i1,i2,i3,i
real(8) v1(3),v2(3),v3(3),t1,t2
complex(8) zt1,zt2
! automatic arrays
real(8) rflm(lmmaxvr)
complex(8) zflm(lmmaxvr)
! external functions
real(8) rfirvec
external rfirvec
! prefactor
zt1=1.d0/(dble(nphcell)*deltaph)
! multiply by i for sin-like displacement
if (iph.eq.1) zt1=zt1*zi
!------------------------------!
!     muffin-tin potential     !
!------------------------------!
ias=0
jas=0
do is=1,nspecies
  ja=0
  do ia=1,natoms0(is)
    ias=ias+1
    do i=1,nphcell
      ja=ja+1
      jas=jas+1
! important: the muffin-tin potential should have an *explicit* phase exp(iq.r)
      t1=-dot_product(vqc(:,iq),atposc(:,ja,is))
      zt2=zt1*cmplx(cos(t1),sin(t1),8)
! loop over radial points
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
! compute the difference between the perturbed and unperturbed potentials
        rflm(:)=veffmt(:,ir,jas)-veffmtp(:,ir,jas)
! convert real potential to a complex spherical harmonic expansion
        call rtozflm(lmaxvr,rflm,zflm)
! add to total
        dveffmt(:,irc,ias)=dveffmt(:,irc,ias)+zt2*zflm(:)
! end loop over radial points
      end do
    end do
! end loop over atoms and species
  end do
end do
!--------------------------------!
!     interstitial potential     !
!--------------------------------!
ir=0
do i3=0,ngrid0(3)-1
  v1(3)=dble(i3)/dble(ngrid0(3))
  do i2=0,ngrid0(2)-1
    v1(2)=dble(i2)/dble(ngrid0(2))
    do i1=0,ngrid0(1)-1
      v1(1)=dble(i1)/dble(ngrid0(1))
      ir=ir+1
      call r3mv(avec0,v1,v2)
      do i=1,nphcell
        v3(:)=v2(:)+vphcell(:,i)
        t1=-dot_product(vqc(:,iq),v3(:))
        zt2=zt1*cmplx(cos(t1),sin(t1),8)
        t1=rfirvec(ngrid,ainv,v3,veffir)
        t2=rfirvec(ngrid,ainv,v3,veffirp)
        dveffir(ir)=dveffir(ir)+zt2*(t1-t2)
      end do
    end do
  end do
end do
return
end subroutine

