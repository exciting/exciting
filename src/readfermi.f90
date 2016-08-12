! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: readfermi
! !INTERFACE:
!
!
subroutine readfermi
! !USES:
  use mod_misc, only: filext
  use mod_eigenvalue_occupancy, only: efermi
! !DESCRIPTION:
!   Reads the Fermi energy from the file {\tt EFERMI.OUT}.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
  implicit none

  open(50, File='EFERMI'//trim(filext), Action='READ',&
   & Form='FORMATTED', Status='OLD')
  read(50,*) efermi
  close(50)

  return
end subroutine
!EOC
