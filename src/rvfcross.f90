
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rvfcross
! !INTERFACE:
subroutine rvfcross(rvfmt1,rvfmt2,rvfir1,rvfir2,rvfmt3,rvfir3)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   rvfmt1 : first input muffin-tin field (in,real(lmmaxvr,nrmtmax,natmtot,3))
!   rvfmt2 : second input muffin-tin field (in,real(lmmaxvr,nrmtmax,natmtot,3))
!   rvfir1 : first input interstitial field (in,real(ngrtot,3))
!   rvfir2 : second input interstitial field (in,real(ngrtot,3))
!   rvfmt3 : output muffin-tin field (out,real(lmmaxvr,nrmtmax,natmtot,3))
!   rvfir3 : output interstitial field (out,real(ngrtot,3))
! !DESCRIPTION:
!   Given two real vector fields, ${\bf f}_1$ and ${\bf f}_2$, defined over the
!   entire unit cell, this routine computes the local cross product
!   $$ {\bf f}_3({\bf r})\equiv{\bf f}_1({\bf r})\times{\bf f}_2({\bf r}). $$
!
! !REVISION HISTORY:
!   Created February 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rvfmt1(lmmaxvr,nrmtmax,natmtot,3)
real(8), intent(in) :: rvfmt2(lmmaxvr,nrmtmax,natmtot,3)
real(8), intent(in) :: rvfir1(ngrtot,3)
real(8), intent(in) :: rvfir2(ngrtot,3)
real(8), intent(out) :: rvfmt3(lmmaxvr,nrmtmax,natmtot,3)
real(8), intent(out) :: rvfir3(ngrtot,3)
! local variables
integer is,ia,ias,ir,itp,i
real(8) v1(3),v2(3),v3(3)
! automatic arrays
real(8) rftp1(lmmaxvr,3),rftp2(lmmaxvr,3)
! muffin-tin region
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      do i=1,3
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
         rvfmt1(:,ir,ias,i),1,0.d0,rftp1(:,i),1)
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
         rvfmt2(:,ir,ias,i),1,0.d0,rftp2(:,i),1)
      end do
      do itp=1,lmmaxvr
        v1(:)=rftp1(itp,:)
        v2(:)=rftp2(itp,:)
        call r3cross(v1,v2,v3)
        rftp1(itp,:)=v3(:)
      end do
      do i=1,3
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp1(:,i),1,0.d0, &
         rvfmt3(:,ir,ias,i),1)
      end do
    end do
  end do
end do
! interstitial region
do ir=1,ngrtot
  v1(:)=rvfir1(ir,:)
  v2(:)=rvfir2(ir,:)
  call r3cross(v1,v2,v3)
  rvfir3(ir,:)=v3(:)
end do
return
end subroutine
!EOC

