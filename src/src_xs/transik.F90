!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
!
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Logical Function transik (ik, trans)
      Use modmain
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik, trans (3, ndftrans)
      transik = .False.
  ! quick return
      If (trans(1, 1) .Eq. 0) Then
         transik = .True.
         Return
      End If
      If (any(trans(1, :) .Eq. ik)) transik = .True.
End Function transik
