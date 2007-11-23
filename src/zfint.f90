
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function zfint(zfmt,zfir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: zfmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) :: zfir(ngrtot)
! local variables
integer is,ia,ias,ir,irc
real(8) t1,t2
complex(8) zsum
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax)
zsum=0.d0
! interstitial contribution
do ir=1,ngrtot
  zsum=zsum+cfunir(ir)*zfir(ir)
end do
zsum=zsum*omega/dble(ngrtot)
! muffin-tin contribution
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do irc=1,nrcmt(is)
      t1=rcmt(irc,is)**2
      fr1(irc)=dble(zfmt(1,irc,ias))*t1
      fr2(irc)=aimag(zfmt(1,irc,ias))*t1
    end do
    call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
    t1=gr(nrcmt(is))
    call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
    t2=gr(nrcmt(is))
    zsum=zsum+fourpi*y00*cmplx(t1,t2,8)
  end do
end do
zfint=zsum
return
end function

