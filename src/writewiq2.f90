!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: writewiq2
! !INTERFACE:
!
!
Subroutine writewiq2
! !USES:
      Use modmain
! !DESCRIPTION:
!   Outputs the integrals of $1/q^2$ in the small parallelepiped around each
!   $q$-point to the file {\tt WIQ2.OUT}. Note that the integrals are calculated
!   after the $q$-point has been mapped to the first Brillouin zone. See routine
!   genwiq2.
!
! !REVISION HISTORY:
!   Created June 2005 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: iq, i1, i2, i3
      Real (8) :: v0 (3), v1 (3), v2 (3), t1, t2
      Open (50, File='WIQ2'//trim(filext), Action='WRITE', Form='FORMAT&
     &TED')
      Write (50, '(I6, " : nqpt; q-point, vql, wiq2 below")') nqpt
      Do iq = 1, nqpt
! map the q-vector into the first Brillouin zone
         t1 = 1.d5
         v0 (:) = 0.d0
         Do i1 = - 1, 1
            Do i2 = - 1, 1
               Do i3 = - 1, 1
                  v1 (:) = vqc (:, iq) + dble (i1) * bvec (:, 1) + dble &
                 & (i2) * bvec (:, 2) + dble (i3) * bvec (:, 3)
                  t2 = v1 (1) ** 2 + v1 (2) ** 2 + v1 (3) ** 2
                  If (t2 .Lt. t1) Then
                     t1 = t2
                     v0 (1) = vql (1, iq) + dble (i1)
                     v0 (2) = vql (2, iq) + dble (i2)
                     v0 (3) = vql (3, iq) + dble (i3)
                     v2 (:) = v1 (:)
                  End If
               End Do
            End Do
         End Do
      End Do
      Close (50)
      Return
End Subroutine
!EOC
