!
!
! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine mapxsparameters
      Use modmain
      Use modinput
      Use modxs
      Implicit None
      input%groundstate%nosym = input%xs%nosym
      input%groundstate%ngridk (:) = input%xs%ngridk(:)
      input%groundstate%reducek = input%xs%reducek
      input%groundstate%vkloff (:) = input%xs%vkloff(:)
      ngridq (:) = input%xs%ngridq(:)
      If (associated(input%phonons)) input%phonons%reduceq = &
     & input%xs%reduceq
      input%groundstate%rgkmax = input%xs%rgkmax
      input%groundstate%swidth = input%xs%swidth
      input%groundstate%lmaxapw = input%xs%lmaxapw
      input%groundstate%lmaxmat = input%xs%lmaxmat
      input%groundstate%nempty = input%xs%nempty
      If (associated(input%groundstate%spin)) Then
         input%groundstate%spin%bfieldc (:) = input%xs%bfieldc(:)
      End If
      input%groundstate%maxscl = input%xs%maxscl
End Subroutine
