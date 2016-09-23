!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findigp0 (ngp, gpc, igp0)
      Implicit None
! arguments
      Integer, Intent (In) :: ngp
      Real (8), Intent (In) :: gpc (ngp)
      Integer, Intent (Out) :: igp0
! local variables
      Integer :: igp
      Real (8), Parameter :: eps = 1.d-14
      Real (8) :: t1
      igp0 = 1
      t1 = gpc (igp0) + eps
      Do igp = 2, ngp
         If (gpc(igp) .Lt. t1) Then
            igp0 = igp
            t1 = gpc (igp) + eps
         End If
      End Do
      Return
End Subroutine
