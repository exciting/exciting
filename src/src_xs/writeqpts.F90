

! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeqpts
! !INTERFACE:


subroutine writeqpts
! !USES:
use modinput
use modmain
use modxs
use m_getunit
use m_genfilname
! !DESCRIPTION:
!   Writes the ${\bf q}$-points in lattice coordinates, weights and number of
!   ${\bf G+q}$-vectors to the file {\tt QPOINTS.OUT}. Based on the routine 
!   {\tt writekpts}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
implicit none
! local variables
integer::iq, un
character(256) :: filnam
call getunit(un)
call genfilname(basename='QPOINTS', appfilext=.true., filnam=filnam)
open(un, file=trim(filnam), action='WRITE', form='FORMATTED')
write(un, '(I6, " : nqpt; q-point, vql, vqc, wqpt, ngq below")') nqpt
do iq=1, nqpt
  write(un, '(I6, 6G18.10, I8)') iq, vql(:, iq), vqc(:, iq), ngq(iq)
end do
close(un)
! write out reduced q-point set for screened Coulomb interaction
if (task.eq.440) then
   call genfilname(basename='QPOINTSR', appfilext=.true., filnam=filnam)
   open(un, file=trim(filnam), action='WRITE', form='FORMATTED', &
	status='replace')
   write(un, '(I6, " : nqptr; q-point, vqlr, vqcr, wqptr below")') nqptr
   do iq=1, nqptr
      write(un, '(I6, 7G18.10)') iq, vqlr(:, iq), vqcr(:, iq), wqptr(iq)
   end do
   close(un)
end if
end subroutine writeqpts
!EOC
