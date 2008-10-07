
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writekpts
! !INTERFACE:
subroutine writekpts
! !USES:
use modmain
! !DESCRIPTION:
!   Writes the $k$-points in lattice coordinates, weights and number of
!   ${\bf G+k}$-vectors to the file {\tt KPOINTS.OUT}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik
open(50,file='KPOINTS'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'(I6," : nkpt; k-point, vkl, wkpt, nmat below")') nkpt
do ik=1,nkpt
  write(50,'(I6,4G18.10,2I8)') ik,vkl(:,ik),wkpt(ik),nmat(:,ik)
end do
close(50)
return
end subroutine
!EOC

