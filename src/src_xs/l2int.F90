!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Integer Function l2int (l)
      Implicit None
      Logical, Intent (In) :: l
      If (l) Then
         l2int = 1
      Else
         l2int = 0
      End If
End Function l2int
