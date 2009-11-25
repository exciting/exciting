!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gentetlinkp
! !INTERFACE:
!
!
Subroutine gentetlinkp (vpl, tqw)
  ! !USES:
      Use modinput
      Use modmain
      Use modtetra
! !DESCRIPTION:
!  Interface to the {\tt gentetlink} routine.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vpl (3)
      Integer, Intent (In) :: tqw
  ! call to interface routine
      Call gentetlink (vpl, tqw, input%structure%epslat, bvec, &
     & input%groundstate%ngridk, input%groundstate%vkloff, nkpt, &
     & nkptnr, vklnr, ikmapnr)
End Subroutine gentetlinkp
!EOC
