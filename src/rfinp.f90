
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rfinp
! !INTERFACE:
real(8) function rfinp(lrstp,rfmt1,rfmt2,rfir1,rfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rfmt1 : first function in real spherical harmonics for all muffin-tins
!           (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfmt2 : second function in real spherical harmonics for all muffin-tins
!           (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfir1 : first real interstitial function in real-space (in,real(ngrtot))
!   rfir2 : second real interstitial function in real-space (in,real(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two real fuctions over the entire unit cell.
!   The input muffin-tin functions should have angular momentum cut-off
!   {\tt lmaxvr}. In the intersitial region, the integrand is multiplied with
!   the characteristic function, $\tilde{\Theta}({\bf r})$, to remove the
!   contribution from the muffin-tin. See routines {\tt rfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(in) :: rfmt1(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: rfmt2(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: rfir1(ngrtot)
real(8), intent(in) :: rfir2(ngrtot)
! local variables
integer is,ia,ias,ir
real(8) sum
! external functions
real(8) rfmtinp
external rfmtinp
sum=0.d0
! interstitial contribution
do ir=1,ngrtot
  sum=sum+rfir1(ir)*rfir2(ir)*cfunir(ir)
end do
sum=sum*omega/dble(ngrtot)
! muffin-tin contribution
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    sum=sum+rfmtinp(lrstp,lmaxvr,nrmt(is),spr(:,is),lmmaxvr,rfmt1(:,:,ias), &
     rfmt2(:,:,ias))
  end do
end do
rfinp=sum
return
end function
!EOC

