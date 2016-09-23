!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine backup1
      Use modmain
      Use modinput
      Use modxs
      Implicit None
      nempty_b = input%groundstate%nempty
      rgkmax_b = input%groundstate%rgkmax
      reducek_b = input%groundstate%reducek
      ngridk_b (:) = input%groundstate%ngridk(:)
      vkloff_b (:) = input%groundstate%vkloff(:)
      emattype_b = input%xs%emattype
End Subroutine backup1
