
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeqpts
! !INTERFACE:
subroutine writeqpts
! !USES:
use modmain
use modtddft
use m_getunit
! !DESCRIPTION:
!   Writes the ${\bf q}$-points in lattice coordinates, weights and number of
!   ${\bf G+q}$-vectors to the file {\tt QPOINTS.OUT}. Based on the routine 
!   "writekpts.f90"
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
implicit none
! local variables
integer iq
call getunit(unit1)
open(unit1,file='QPOINTS'//trim(filext),action='WRITE',form='FORMATTED')
write(unit1,'(I6," : nqpt; q-point, vql, vqc, wqpt, ngq below")') nqpt
do iq=1,nqpt
  write(unit1,'(I6,7G18.10,2I8)') iq,vql(:,iq),vqc(:,iq),wqpt(iq), ngq(iq)
end do
close(unit1)
end subroutine writeqpts
!EOC

