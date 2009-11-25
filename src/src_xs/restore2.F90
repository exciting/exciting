!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine restore2
      Use modmain
      Use modinput
      Use modxs
      Implicit None
      ngridq (:) = ngridq_b (:)
      If (associated(input%phonons)) Then
         input%phonons%reduceq = reduceq_b
      End If
End Subroutine restore2
