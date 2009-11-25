!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: r3frac
! !INTERFACE:
!
!
Subroutine r3frac (eps, v, iv)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
!   iv  : integer parts of v (out,integer(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!   The integer components of {\tt v} are returned in the variable {\tt iv}.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: eps
      Real (8), Intent (Inout) :: v (3)
      Integer, Intent (Out) :: iv (3)
! local variables
      Integer :: i
      Do i = 1, 3
         iv (i) = Int (v(i))
         v (i) = v (i) - dble (iv(i))
         If (v(i) .Lt. 0.d0) Then
            v (i) = v (i) + 1.d0
            iv (i) = iv (i) - 1
         End If
         If (1.d0-v(i) .Lt. eps) Then
            v (i) = 0.d0
            iv (i) = iv (i) + 1
         End If
         If (v(i) .Lt. eps) Then
            v (i) = 0.d0
         End If
      End Do
      Return
End Subroutine
!EOC
