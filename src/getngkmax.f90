
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getngkmax
! !INTERFACE:
subroutine getngkmax
! !USES:
use modmain
! !DESCRIPTION:
!   Determines the largest number of ${\bf G+k}$-vectors with length less than
!   {\tt gkmax} over all the $k$-points and stores it in the global variable
!   {\tt ngkmax}. This variable is used for allocating arrays.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ispn,ik,i,ig
real(8) v1(3),v2(3),t1,t2
t1=gkmax**2
ngkmax=0
do ispn=1,nspnfv
  do ik=1,nkpt
    if (spinsprl) then
! spin-spiral case
      if (ispn.eq.1) then
        v1(:)=vkc(:,ik)+0.5d0*vqcss(:)
      else
        v1(:)=vkc(:,ik)-0.5d0*vqcss(:)
      end if
    else
      v1(:)=vkc(:,ik)
    end if
    i=0
    do ig=1,ngvec
      v2(:)=vgc(:,ig)+v1(:)
      t2=v2(1)**2+v2(2)**2+v2(3)**2
      if (t2.lt.t1) i=i+1
    end do
    ngkmax=max(ngkmax,i)
  end do
end do
return
end subroutine
!EOC
