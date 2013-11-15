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
      Integer :: is, ia, fsn, i, ias

      fsn = input%groundstate%spin%fixspinnumber

      if (fsn .eq. 0) return

      i = 0
      Write (fnum,*)


      if ( (fsn .eq. 1) .or. (fsn .eq. 3) ) then
         write (fnum, '(" FSM global effective field", T35, ":", 3F15.8)') bfsmc(1:ndmag)
      end if
      if ( (fsn .eq. 2) .or. (fsn .eq. 3) ) then
         write (fnum, '(" FSM effective field in muffin-tin:")')
         do is = 1, nspecies
            do ia = 1, natoms (is)
               ias = idxas (ia, is)
               i = i+1
               write(fnum,'("     atom ",I5,4x,A2,T35,":",3F15.8)') &
              &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), bfsmcmt(1:ndmag,ia,is)
            end do
         end do
      end if

      Return
End Subroutine
