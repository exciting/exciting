
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: charge
! !INTERFACE:
subroutine charge
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total charges by integrating the
!   density.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir
real(8) sum,t1
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
! find the muffin-tin charges
chgmttot=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      fr(ir)=rhomt(1,ir,ias)*spr(ir,is)**2
    end do
    call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
    chgmt(ias)=fourpi*y00*gr(nrmt(is))
    chgmttot=chgmttot+chgmt(ias)
  end do
end do
! find the interstitial charge
sum=0.d0
do ir=1,ngrtot
  sum=sum+rhoir(ir)*cfunir(ir)
end do
chgir=sum*omega/dble(ngrtot)
! total calculated charge
chgcalc=chgmttot+chgir
t1=chgtot/chgcalc
if (abs(t1-1.d0).gt.epschg) then
  write(*,*)
  write(*,'("Warning(charge): total charge density incorrect for s.c. &
   &loop ",I5)') iscl
  write(*,'(" Calculated : ",G18.10)') chgcalc
  write(*,'(" Required   : ",G18.10)') chgtot
end if
return
end subroutine
!EOC
