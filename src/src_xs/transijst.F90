!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Logical Function transijst (ik, ist, jst, trans)
      Use modmain
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik, ist, jst, trans (3, ndftrans)
  ! local variables
      Integer :: l, ikt, it, jt, ir, jr
      transijst = .False.
      Do l = 1, ndftrans
         ikt = trans (1, l)
         it = trans (2, l)
         jt = trans (3, l)
         If ((ikt .Eq. 0) .Or. (ikt .Eq. ik)) Then
        ! it-jt, it-x, x-jt and x-x combinations
            If (((it .Eq. 0) .Or. (it .Eq. ist)) .And. ((jt .Eq. 0) &
           & .Or. (jt .Eq. jst))) Then
               transijst = .True.
               Return
            End If
        ! ranges between |it| and |jt| (for negative it and jt)
            ir = - it
            jr = - jt
            If (((ir .Gt. 0) .And. (jr .Gt. 0)) .And. (ir .Lt. jr)) &
           & Then
               If ((ist .Ge. ir) .And. (jst .Le. jr)) transijst = &
              & .True.
               Return
            End If
         End If
      End Do
End Function transijst
