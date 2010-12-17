
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initsolver
      use modinput
      Use inputdom
      Implicit None
      If ( .Not. (associated(input%groundstate%solver))) Then
        ! set the default values if solver element not present
         input%groundstate%solver => getstructsolver (emptynode)
      End If
end subroutine
