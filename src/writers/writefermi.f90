!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writefermi
! !INTERFACE:
!
!
Subroutine writefermi
! !USES:
      Use modmain
! !DESCRIPTION:
!   Writes the Fermi energy to the file {\tt EFERMI.OUT}.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
      Implicit None
      Open (50, File='EFERMI'//trim(filext), Action='WRITE', Form='FORM&
     &ATTED')
      Write (50, '(G18.10)') efermi
      Close (50)
      Return
End Subroutine
!EOC
