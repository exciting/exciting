
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine tasklauncher
      Use modinput
      Use modmain, Only: task
      Use inputdom
!
      Implicit None
!
      If (associated(input%groundstate)) Then
         If ( .Not. (associated(input%groundstate%solver))) Then
  		! set the default values if solver element not present
            input%groundstate%solver => getstructsolver (emptynode)
         End If
         If (input%groundstate%do .Eq. "fromscratch") Then
            If (associated(input%structureoptimization)) Then
               task = 2
            Else
               task = 0
            End If
         Else
            If (associated(input%structureoptimization)) Then
               task = 3
            Else
               task = 1
            End If
         End If
         If (input%groundstate%do .Ne. "skip") Call gndstate
      End If
!
!
      If (associated(input%properties)) Then
         Call propertylauncher
      End If
!
      If (associated(input%phonons)) Then
         call phononstasklauncher()
      End If
!
      If (associated(input%xs)) Then
         Call xstasklauncher ()
      End If
End Subroutine
