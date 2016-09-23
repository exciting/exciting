
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeqpts
! !INTERFACE:
subroutine gw_writeqpts
! !USES:
use modmain
use modgw
! !DESCRIPTION:
!   Writes the ${\bf q}$-points in lattice coordinates, weights and number of
!   ${\bf G+q}$-vectors to the file {\tt QPOINTS.OUT}.
!
! !REVISION HISTORY:
!   Created May 2006 (RGA)
!EOP
!BOC
implicit none
! local variables
integer iq
open(50,file='QPOINTS'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'(I6," : nqpt; q-point, vql, wqpt, ngq below")') nqpt
do iq=1,nqpt
  write(50,'(I6,4G18.10,2I8)') iq,vql(:,iq),wkpt(iq),ngq(iq)
end do
close(50)
return
end subroutine
!EOC

