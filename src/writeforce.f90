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
      Write (fnum,'("Total atomic forces:")')
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            write(fnum,'(" atom ",I4,4x,A2,T30,": ",3F14.8)') &
           &  ia, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
           &  forcetot(:,ias)
            write(fnum, '(T30,"Hellmann-Feynman",4x,3F14.8)') forcehf(:,ias)
            write(fnum, '(T30,"core correction",4x,3F14.8)') forcecr(:,ias)
            write (fnum, '("   IBS", T30, ": ", 3F14.8)') forceibs (:,ias)
            t1 = Sqrt (forcetot(1, ias)**2+forcetot(2,ias)**2+forcetot(3, ias)**2)
            write(fnum,'(" total magnitude",T30,": ",F14.8)') t1
         End Do
      End Do
      Return
End Subroutine
