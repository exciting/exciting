!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine chkptchk
      Use modxs
      Implicit None
  ! local variables
      Logical :: exis
      Inquire (File=trim(fnresume), Exist=exis)
      If (exis) Then
         Write (*,*)
         Write (*, '("Error(chkptchk): stale checkpoint file found.")')
         Write (*, '(" Either your previous calculation crashed or anot&
        &her")')
         Write (*, '(" instance of exciting is already running")')
         Write (*,*)
         Call terminate
      End If
End Subroutine chkptchk
