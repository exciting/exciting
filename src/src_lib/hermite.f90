!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: hermite
! !INTERFACE:
Real (8) Function hermite (n, x)
! !INPUT/OUTPUT PARAMETERS:
!   n : order of Hermite polynomial (in,integer)
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the $n$th Hermite polynomial. The recurrence relation
!   $$ H_i(x)=2xH_{i-1}(x)-2nH_{i-2}(x), $$
!   with $H_0=1$ and $H_1=2x$, is used. This procedure is numerically stable
!   and accurate to near machine precision for $n\le 20$.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: x
! local variables
      Integer :: i
      Real (8) :: h1, h2, ht
      If (n .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(hermite): n < 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (n .Gt. 20) Then
         Write (*,*)
         Write (*, '("Error(hermite): n out of range : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (Abs(x) .Gt. 1.d15) Then
         Write (*,*)
         Write (*, '("Error(hermite): x out of range : ", G18.10)') x
         Write (*,*)
         Stop
      End If
      If (n .Eq. 0) Then
         hermite = 1.d0
         Return
      End If
      If (n .Eq. 1) Then
         hermite = 2.d0 * x
         Return
      End If
      h1 = 2.d0 * x
      h2 = 1.d0
      Do i = 2, n
         ht = 2.d0 * x * h1 - 2.d0 * dble (i-1) * h2
         h2 = h1
         h1 = ht
      End Do
      hermite = h1
      Return
End Function
!EOC
