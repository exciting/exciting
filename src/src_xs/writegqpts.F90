!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writegqpts
      Implicit None
Contains
!
!BOP
! !ROUTINE: writegqpts
! !INTERFACE:
!
!
      Subroutine writegqpts (iq, filex)
! !USES:
         Use modmain
         Use modxs
         Use m_getunit
! !DESCRIPTION:
!   Writes the ${\bf G+q}$-points in lattice coordinates, Cartesian
!   coordinates, and lengths of ${\bf G+q}$-vectors to the file
!   {\tt GQPOINTS.OUT}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq
         Character (*), Intent (In) :: filex
    ! local variables
         Integer :: igq
         Call getunit (unit1)
         Open (unit1, File='GQPOINTS'//trim(filex), Action='WRITE', &
        & Form='FORMATTED')
         Write (unit1, '(I6, " : ngq; G+q-point, vql, vqc, wqpt, ngq be&
        &low")') ngq (iq)
         Do igq = 1, ngq (iq)
            Write (unit1, '(I6, 7G18.10)') igq, vgql (:, igq, iq), vgqc &
           & (:, igq, iq), gqc (igq, iq)
         End Do
         Close (unit1)
      End Subroutine writegqpts
!EOC
!
End Module m_writegqpts
