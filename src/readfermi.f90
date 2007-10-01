
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readfermi
! !INTERFACE:
subroutine readfermi
! !USES:
use modmain
! !DESCRIPTION:
!   Reads the Fermi energy from the file {\tt EFERMI.OUT}.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
open(50,file='EFERMI'//trim(filext),action='READ',form='FORMATTED',status='OLD')
read(50,*) efermi
close(50)
return
end subroutine
!EOC

