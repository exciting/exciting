!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writekpts
! !INTERFACE:
!
!
Subroutine writekpts
! !USES:
      Use modmain
! !DESCRIPTION:
!   Writes the $k$-points in lattice coordinates, weights and number of
!   ${\bf G+k}$-vectors to the file {\tt KPOINTS.OUT}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik
      Open (50, File='KPOINTS'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Write (50, '(I6, " : nkpt; k-point, vkl, wkpt, nmat below")') &
     & nkpt
      Do ik = 1, nkpt
         Write (50, '(I6, 4G18.10, 2I8)') ik, vkl (:, ik), wkpt (ik), &
        & nmat (:, ik)
      End Do
      Close (50)
      Return
End Subroutine
!EOC
