!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writeqpts
! !INTERFACE:
!
!
Subroutine writeqpts
! !USES:
      Use modinput
      Use modmain
      Use modxs
      Use m_getunit
      Use m_genfilname
! !DESCRIPTION:
!   Writes the ${\bf q}$-points in lattice coordinates, weights and number of
!   ${\bf G+q}$-vectors to the file {\tt QPOINTS.OUT}. Based on the routine
!   {\tt writekpts}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: iq, un
      Character (256) :: filnam
      Call getunit (un)
      Call genfilname (basename='QPOINTS', appfilext=.True., &
     & filnam=filnam)
      Open (un, File=trim(filnam), Action='WRITE', Form='FORMATTED')
      Write (un, '(I6, " : nqpt; q-point, vql, vqc, wqpt, ngq below")') &
     & nqpt
      Do iq = 1, nqpt
         Write (un, '(I6, 6G18.10, I8)') iq, vql (:, iq), vqc (:, iq), &
        & ngq (iq)
      End Do
      Close (un)
! write out reduced q-point set for screened Coulomb interaction
      If (task .Eq. 440) Then
         Call genfilname (basename='QPOINTSR', appfilext=.True., &
        & filnam=filnam)
         Open (un, File=trim(filnam), Action='WRITE', Form='FORMATTED', &
        & Status='replace')
         Write (un, '(I6, " : nqptr; q-point, vqlr, vqcr, wqptr below")&
        &') nqptr
         Do iq = 1, nqptr
            Write (un, '(I6, 7G18.10)') iq, vqlr (:, iq), vqcr (:, iq), &
           & wqptr (iq)
         End Do
         Close (un)
      End If
End Subroutine writeqpts
!EOC
