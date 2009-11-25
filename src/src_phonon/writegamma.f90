!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writegamma (gq)
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: gq (3*natmtot, nqpt)
! local variables
      Integer :: iq, i
      Open (50, File='GAMMAQ.OUT', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '(I4, " : total number of atoms")') natmtot
      Write (50, '(I6, " : number of q-points")') nqpt
      Write (50,*)
      Do iq = 1, nqpt
         Write (50, '(I6, " : q-point")') iq
         Write (50, '(3G18.10, " : q-vector (lattice coordinates)")') &
        & vql (:, iq)
         Write (50, '(3G18.10, " : q-vector (Cartesian coordinates)")') &
        & vqc (:, iq)
         Do i = 1, 3 * natmtot
            Write (50, '(I4, G18.10)') i, gq (i, iq)
         End Do
         Write (50,*)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writegamma):")')
      Write (*, '(" wrote phonon linewidths for all q-points to GAMMAQ.&
     &OUT")')
      Write (*,*)
      Return
End Subroutine
