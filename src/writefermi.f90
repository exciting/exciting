
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writefermi
! !INTERFACE:
subroutine writefermi
! !USES:
use modmain
! !DESCRIPTION:
!   Writes the Fermi energy to the file {\tt EFERMI.OUT}.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
open(50,file='EFERMI'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'(G18.10)') efermi
close(50)
return
end subroutine
!EOC

