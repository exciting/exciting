!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writelinen
! !INTERFACE:
!
!
Subroutine writelinen
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Writes the linearisation energies for all APW and local-orbital functions to
!   the file {\tt LINENGY.OUT}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, l, io, ilo
      Open (50, File='LINENGY'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Write (50,*)
            Write (50, '("Species : ", I4, " (", A, "), atom : ", I4)') &
           & is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol), &
           & ia
            Write (50, '(" APW functions :")')
            Do l = 0, input%groundstate%lmaxapw
               Do io = 1, apword (l, is)
                  Write (50, '("  l = ", I2, ", order = ", I2, " : ", G&
                 &18.10)') l, io, apwe (io, l, ias)
               End Do
            End Do
            Write (50, '(" local-orbital functions :")')
            Do ilo = 1, nlorb (is)
               Do io = 1, lorbord (ilo, is)
                  Write (50, '("  l.o. = ", I2, ", l = ", I2, ", order &
                 &= ", I2, " : ", G18.10)') ilo, lorbl (ilo, is), io, &
                 & lorbe (io, ilo, ias)
               End Do
            End Do
         End Do
      End Do
      Close (50)
      Return
End Subroutine
!EOC
