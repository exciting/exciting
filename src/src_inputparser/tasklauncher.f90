
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine tasklauncher
      Use modinput
      Implicit None

      If (associated(input%groundstate)) &
        call groundstatetasklauncher()

      If (associated(input%properties)) &
        Call propertylauncher

      If (associated(input%phonons)) &
        call phononstasklauncher()

      If (associated(input%xs)) &
        Call xstasklauncher ()
End Subroutine
