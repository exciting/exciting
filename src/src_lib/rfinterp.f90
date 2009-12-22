!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rfinterp
! !INTERFACE:
!
!
Subroutine rfinterp (ni, xi, ldi, fi, no, xo, ldo, fo)
! !INPUT/OUTPUT PARAMETERS:
!   ni  : number of input points (in,integer)
!   xi  : input abscissa array (in,real(ni))
!   ldi : leading dimension (in,integer)
!   fi  : input data array (in,real(ldi,ni)
!   no  : number of output points (in,integer)
!   xo  : output abscissa array (in,real(ni))
!   ldo : leading dimension (in,integer)
!   fo  : output interpolated function (out,real(ldo,no))
! !DESCRIPTION:
!   Given a function defined on a set of input points, this routine uses a
!   clamped cubic spline to interpolate the function on a different set of
!   points. See routine {\tt spline}.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!EOP
!BOC
      Implicit None
      Integer, Intent (In) :: ni
      Real (8), Intent (In) :: xi (ni)
      Integer, Intent (In) :: ldi
      Real (8), Intent (In) :: fi (ldi, ni)
      Integer, Intent (In) :: no
      Real (8), Intent (In) :: xo (no)
      Integer, Intent (In) :: ldo
      Real (8), Intent (Out) :: fo (ldo, no)
! local variables
      Integer :: i, j, k, l
      Real (8) :: dx, t
! automatic arrays
      Real (8) :: cf (3, ni)
      If (ni .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rfinterp): invalid number of input points :&
        & ", I8)') ni
         Write (*,*)
         Stop
      End If
      If (no .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rfinterp): invalid number of output points : ", I8)') no
         Write (*,*)
         Stop
      End If
      If (ni .Eq. 1) Then
         fo (1, :) = fi (1, 1)
         Return
      End If
      Call spline (ni, xi, ldi, fi, cf)
! evaluate spline at output points
      i = 1
      Do l = 1, no
         t = xo (l)
         If (i .Ge. ni) i = 1
         If (t .Lt. xi(i)) Go To 10
         If (t .Le. xi(i+1)) Go To 30
! binary search
10       Continue
         i = 1
         j = ni + 1
20       Continue
         k = (i+j) / 2
         If (t .Lt. xi(k)) Then
            j = k
         Else
            i = k
         End If
         If (j .Gt. i+1) Go To 20
30       Continue
         dx = t - xi (i)
         fo (1, l) = fi (1, i) + dx * (cf(1, i)+dx*(cf(2, i)+dx*cf(3, &
        & i)))
      End Do
      Return
End Subroutine
!EOC
