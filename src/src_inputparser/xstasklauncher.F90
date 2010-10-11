
! Copyright (C) 2009-2010 S. Sagmeister, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine xstasklauncher
      Use modinput
      Use modmain, Only: task
      Use inputdom
!
      If ( .Not. (associated(input%xs%tddft))) Then
  ! set the default values if tddft element not present
         input%xs%tddft => getstructtddft (emptynode)
      End If
!
      If ( .Not. (associated(input%xs%screening))) Then
  ! set the default values if screening element not present
         input%xs%screening => getstructscreening (emptynode)
      End If
!
      If ( .Not. (associated(input%xs%BSE))) Then
  ! set the default values if BSE element not present
         input%xs%BSE => getstructBSE (emptynode)
      End If
!
      If ( .Not. (associated(input%xs%tetra))) Then
  ! set the default values if tetra element not present
         input%xs%tetra => getstructtetra (emptynode)
      End If
!
      Call backup0
      Call backup1
      Call backup2
!
      If (associated(input%xs%plan)) Then
         Call xsmain (input%xs%plan)
      Else If (trim(input%xs%xstype) .Eq. "TDDFT") Then
!
         If (input%xs%tddft%do .eq. "fromkernel") Go To 10
!
         task = 301
         Call xsinit
         Call xsgeneigvec
         Call xsfinit
!
         If ((input%xs%tetra%tetradf)) Then
            task = 310
            Call xsinit
            Call tetcalccw
            Call xsfinit
         End If
!
         task = 320
         Call xsinit
         Call writepmatxs
         Call xsfinit
!
         task = 330
         Call xsinit
         Call writeemat
         Call xsfinit
!
         If (input%xs%tddft%fxctypenumber .Eq. 7 .Or. &
        & input%xs%tddft%fxctypenumber .Eq. 8) Then
            task = 401
            Call xsinit
            Call scrgeneigvec
            Call xsfinit
!
            task = 420
            Call xsinit
            Call scrwritepmat
            Call xsfinit
!
            If ((input%xs%tetra%tetradf)) Then
               task = 410
               Call xsinit
               Call scrtetcalccw
               Call xsfinit
            End If
!
            If (input%xs%screening%do .Eq. "fromscratch") Then
               task = 430
               Call xsinit
               Call screen
               Call xsfinit
!
               task = 440
               Call xsinit
               Call scrcoulint
               Call xsfinit
            End If
!
            task = 450
            Call xsinit
            Call kernxc_bse
            Call xsfinit
         End If
!
         task = 340
         Call xsinit
         Call df
         Call xsfinit
!
10       Continue
!
         task = 350
         Call xsinit
         Call idf
         Call xsfinit
!
      Else If (trim(input%xs%xstype) .Eq. "BSE") Then
!
         task = 301
         Call xsinit
         Call xsgeneigvec
         Call xsfinit
!
         task = 320
         Call xsinit
         Call writepmatxs
         Call xsfinit
!
         task = 401
         Call xsinit
         Call scrgeneigvec
         Call xsfinit
!
         If ((input%xs%tetra%tetradf)) Then
            task = 410
            Call xsinit
            Call scrtetcalccw
            Call xsfinit
         End If
!
         task = 420
         Call xsinit
         Call scrwritepmat
         Call xsfinit
!
         If (input%xs%screening%do .Eq. "fromscratch") Then
            task = 430
            Call xsinit
            Call screen
            Call xsfinit
!
            task = 440
            Call xsinit
            Call scrcoulint
            Call xsfinit
         End If
!
         task = 441
         Call xsinit
         Call exccoulint
         Call xsfinit
!
         task = 445
         Call xsinit
         Call BSE
         Call xsfinit
      Else
         Write (*,*) "error xstasklauncher"
         Write (*,*) trim (input%xs%xstype), "no valid xstype"
         Stop
      End If
!
!
End Subroutine
