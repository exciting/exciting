
! Copyright (C) 2009-2010 S. Sagmeister, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine gwtasklauncher
      Use modinput
      Use modmain, Only: task
      Use inputdom

      if (associated(input%gw)) then
         task = 1
         call gw_main
      else
         write (*,*) "error gwtasklauncher"
         stop
      end if

End Subroutine gwtasklauncher
