
! Copyright (C) 2002-2005 S. Sagmeister J. K. Dewhurst, S. Sharma
! and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writegqpts
  implicit none
contains

!BOP
! !ROUTINE: writegqpts
! !INTERFACE:
subroutine writegqpts(iq)
! !USES:
use modmain
use modtddft
use m_getunit
! !DESCRIPTION:
!   Writes the ${\bf G+q}$-points in lattice coordinates, Cartesian 
!   coordinates, and lengths of ${\bf G+q}$-vectors to the file 
!   {\tt QPOINTS.OUT}. Based on the routine {\tt writekpts.f90}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: iq
! local variables
integer :: igq
call getunit(unit1)
open(unit1,file='GQPOINTS'//trim(filext),action='WRITE',form='FORMATTED')
write(unit1,'(I6," : ngq; G+q-point, vql, vqc, wqpt, ngq below")') ngq(iq)
do igq=1,ngq(iq)
   write(unit1,'(I6,7G18.10)') igq, vgql(:,igq,iq), vgqc(:,igq,iq), &
        gqc(igq,iq)
end do
close(unit1)
end subroutine writegqpts
!EOC

end module m_writegqpts
