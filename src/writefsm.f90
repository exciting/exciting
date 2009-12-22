!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writefsm (fnum)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
! local variables
      Integer :: is, ia
      If (input%groundstate%spin%fixspinnumber .Eq. 0) Return
      Write (fnum,*)
      If ((input%groundstate%spin%fixspinnumber .Eq. 1) .Or. &
     & (input%groundstate%spin%fixspinnumber .Eq. 3)) Then
         Write (fnum, '("FSM global effective field", T30, ": ", 3G18.1&
        &0)') bfsmc (1:ndmag)
      End If
      If ((input%groundstate%spin%fixspinnumber .Eq. 2) .Or. &
     & (input%groundstate%spin%fixspinnumber .Eq. 3)) Then
         Write (fnum, '("FSM local muffin-tin effective fields :")')
         Do is = 1, nspecies
            Write (fnum, '(" species : ", I4, " (", A, ")")') is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol)
            Do ia = 1, natoms (is)
               Write (fnum, '("  atom ", I4, T30, ": ", 3G18.10)') ia, &
              & bfsmcmt (1:ndmag, ia, is)
            End Do
         End Do
      End If
      Return
End Subroutine
