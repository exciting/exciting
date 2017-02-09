!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematbdlims (typ, n1, lo1, hi1, n2, lo2, hi2)
      Use modmain
      Use modxs
      Use modxas
      Implicit None
  ! arguments
      Integer, Intent (In) :: typ
      Integer, Intent (Out) :: n1, n2, lo1, hi1, lo2, hi2
	If (task .GE. 440 .AND. task .LE. 446) Then
      Select Case (typ)
      Case (0)
     ! all combinations
         lo1 = 1
         hi1 = nstsv
         lo2 = 1
         hi2 = nstsv
         n1 = nstsv
         n2 = nstsv
      Case (1)
     ! o-u combinations
         lo1 = sta1
         hi1 = sto1
         lo2 = istunocc0+sta2-1
         hi2 = istunocc0+sto2-1
         n1 = sto1- sta1 + 1
         n2 = sto2 -sta2 + 1
      Case (2)
     ! u-o combinations
         lo1 = istunocc0+sta2-1
         hi1 = istunocc0+stO2-1
         lo2 = sta1
         hi2 = sto1
         n1 = sto2-sta2+1
         n2 = sto1-sta1+1
      Case (3)
     ! o-o combinations
         lo1 = sta1
         hi1 = sto1
         lo2 = sta1
         hi2 = sto1
         n1 = sto1-sta1+1
         n2 = sto1-sta1+1
      Case (4)
     ! u-u combinations
         lo1 = istunocc0+sta2-1
         hi1 = istunocc0+sto2-1
         lo2 = istunocc0+sta2-1
         hi2 = istunocc0+sto2-1
         n1 = sto2-sta2+1
         n2 = sto2-sta2+1
      End Select
     Else
      Select Case (typ)
      Case (0)
     ! all combinations
         lo1 = 1
         hi1 = nstsv
         lo2 = 1
         hi2 = nstsv
         n1 = nstsv
         n2 = nstsv
      Case (1)
     ! o-u combinations
         lo1 = 1
         hi1 = istocc0
         lo2 = istunocc0
         hi2 = nstsv
         n1 = istocc0 - 1 + 1
         n2 = nstsv - istunocc0 + 1
      Case (2)
     ! u-o combinations
         lo1 = istunocc0
         hi1 = nstsv
         lo2 = 1
         hi2 = istocc0
         n1 = nstsv - istunocc0 + 1
         n2 = istocc0 - 1 + 1
      Case (3)
     ! o-o combinations
         lo1 = 1
         hi1 = istocc0
         lo2 = 1
         hi2 = istocc0
         n1 = istocc0 - 1 + 1
         n2 = istocc0 - 1 + 1
      Case (4)
     ! u-u combinations
         lo1 = istunocc0
         hi1 = nstsv
         lo2 = istunocc0
         hi2 = nstsv
         n1 = nstsv - istunocc0 + 1
         n2 = nstsv - istunocc0 + 1
      End Select
	End If
End Subroutine ematbdlims
