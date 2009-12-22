!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sphcrd
! !INTERFACE:
!
!
Subroutine sphcrd (v, r, tp)
! !INPUT/OUTPUT PARAMETERS:
!   v  : input vector (in,real(3))
!   r  : length of v (out,real)
!   tp : (theta, phi) coordinates (out,real(2))
! !DESCRIPTION:
!   Returns the spherical coordinates $(r,\theta,\phi)$ of a vector
!   $$ {\bf v}=(r\sin(\theta)\cos(\phi), r\sin(\theta)\sin(\phi),
!    r\cos(\theta)). $$
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: v (3)
      Real (8), Intent (Out) :: r
      Real (8), Intent (Out) :: tp (2)
! local variables
      Real (8), Parameter :: twopi = 6.2831853071795864769d0
      Real (8), Parameter :: eps = 1.d-14
      Real (8) :: t1
      r = Sqrt (v(1)**2+v(2)**2+v(3)**2)
      If (r .Gt. eps) Then
         t1 = v (3) / r
         If (t1 .Gt. 0.d0) Then
            t1 = Min (t1, 1.d0)
         Else
            t1 = Max (t1,-1.d0)
         End If
         tp (1) = Acos (t1)
         If ((Abs(v(1)) .Gt. eps) .Or. (Abs(v(2)) .Gt. eps)) Then
            tp (2) = Atan2 (v(2), v(1))
            If (tp(2) .Lt. 0.d0) tp (2) = tp (2) + twopi
         Else
            tp (2) = 0.d0
         End If
      Else
         tp (1) = 0.d0
         tp (2) = 0.d0
      End If
      Return
End Subroutine
!EOC
