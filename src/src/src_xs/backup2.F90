!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine backup2
      Use modmain
      Use modinput
      Use modxs
      Implicit None
      ngridq_b (:) = ngridq (:)
      If (associated(input%phonons)) Then
         reduceq_b = input%phonons%reduceq
      End If
End Subroutine backup2
