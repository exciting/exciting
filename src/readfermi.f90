!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: readfermi
! !INTERFACE:
!
!
Subroutine readfermi
! !USES:
      Use modmain
! !DESCRIPTION:
!   Reads the Fermi energy from the file {\tt EFERMI.OUT}.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
      Implicit None
      Open (50, File='EFERMI'//trim(filext), Action='READ', Form='FORMA&
     &TTED', Status='OLD')
      Read (50,*) efermi
      Close (50)
      Return
End Subroutine
!EOC
