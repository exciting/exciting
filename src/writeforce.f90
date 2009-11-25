!
!
!
! Copyright (C) 2004-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writeforce (fnum)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
! local variables
      Integer :: is, ia, ias
      Real (8) :: t1
      Write (fnum,*)
      Write (fnum, '("Forces :")')
      Do is = 1, nspecies
         Write (fnum, '(" species : ", I4, " (", A, ")")') is, trim &
        & (input%structure%speciesarray(is)%species%chemicalSymbol)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Write (fnum, '("  atom : ", I4)') ia
            Write (fnum, '("   Hellmann-Feynman", T30, ": ", 3F14.8)') &
           & forcehf (:, ias)
            Write (fnum, '("   core correction", T30, ": ", 3F14.8)') &
           & forcecr (:, ias)
            Write (fnum, '("   IBS", T30, ": ", 3F14.8)') forceibs (:, &
           & ias)
            Write (fnum, '("   total force", T30, ": ", 3F14.8)') &
           & forcetot (:, ias)
            t1 = Sqrt (forcetot(1, ias)**2+forcetot(2, &
           & ias)**2+forcetot(3, ias)**2)
            Write (fnum, '("   total magnitude", T30, ": ", F14.8)') t1
         End Do
      End Do
!call flushifc(fnum)
      Return
End Subroutine
