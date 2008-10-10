
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zrhogp(gpc,jlgpr,ylmgp,sfacgp,zrhomt,zrhoir,zrho0)
use modmain
implicit none
! arguments
real(8), intent(in) :: gpc
real(8), intent(in) :: jlgpr(0:lmaxvr,nrcmtmax,nspecies)
complex(8), intent(in) :: ylmgp(lmmaxvr)
complex(8), intent(in) :: sfacgp(natmtot)
complex(8), intent(in) :: zrhomt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) :: zrhoir(ngrtot)
complex(8), intent(out) :: zrho0
! local variables
integer is,ia,ias
integer l,m,lm,nrc,ir
real(8) t1,t2
complex(8) zsum1,zsum2
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
real(8) gr(nrcmtmax),cf(3,nrcmtmax)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
! (note that the phase exp(ip.r) is implicit)
zrho0=0.d0
do ir=1,ngrtot
  zrho0=zrho0+cfunir(ir)*zrhoir(ir)
end do
zrho0=zrho0/dble(ngrtot)
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
! (note that the phase exp(ip.r) is explicit)
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrc
      zsum1=0.d0
      lm=0
      do l=0,lmaxvr
        lm=lm+1
        zsum2=zrhomt(lm,ir,ias)*ylmgp(lm)
        do m=-l+1,l
          lm=lm+1
          zsum2=zsum2+zrhomt(lm,ir,ias)*ylmgp(lm)
        end do
        zsum1=zsum1+jlgpr(l,ir,is)*conjg(zil(l))*zsum2
      end do
      t1=rcmt(ir,is)**2
      fr1(ir)=dble(zsum1)*t1
      fr2(ir)=aimag(zsum1)*t1
    end do
    call fderiv(-1,nrc,rcmt(:,is),fr1,gr,cf)
    t1=gr(nrc)
    call fderiv(-1,nrc,rcmt(:,is),fr2,gr,cf)
    t2=gr(nrc)
    zrho0=zrho0+(fourpi/omega)*conjg(sfacgp(ias))*cmplx(t1,t2,8)
  end do
end do
return
end subroutine

