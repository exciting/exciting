! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! main routine for the EXCITING code
!
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
         Call phonon
      End If
!
      If (associated(input%xs)) Then
         Call xstasklauncher ()
      End If
End Subroutine
