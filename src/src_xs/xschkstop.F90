!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xschkstop
      Use modxs
      use modmpi
      Use m_filedel
      Implicit None
  ! local variables
      Logical :: exist
      Inquire (File='STOP', Exist=Exist)
      If (exist) Then
         Write (unitout, '("STOP file exists - stopping with message: "&
        &, a)') trim (msg)
         Call filedel ('STOP')
         Call terminate
      End If
End Subroutine xschkstop
