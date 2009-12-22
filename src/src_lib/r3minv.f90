!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: r3minv
! !INTERFACE:
!
!
Subroutine r3minv (a, b)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,real(3,3))
!   b : output matrix (in,real(3,3))
! !DESCRIPTION:
!   Computes the inverse of a real $3\times 3$ matrix.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: a (3, 3)
      Real (8), Intent (Out) :: b (3, 3)
! local variables
      Real (8) :: t1
      t1 = a (1, 2) * a (2, 3) * a (3, 1) - a (1, 3) * a (2, 2) * a (3, &
     & 1) + a (1, 3) * a (2, 1) * a (3, 2) - a (1, 1) * a (2, 3) * a &
     & (3, 2) + a (1, 1) * a (2, 2) * a (3, 3) - a (1, 2) * a (2, 1) * &
     & a (3, 3)
      If (Abs(t1) .Lt. 1.d-40) Then
         Write (*,*)
         Write (*, '("Error(r3minv): singular matrix")')
         Write (*,*)
         Stop
      End If
      t1 = 1.d0 / t1
      b (1, 1) = (a(2, 2)*a(3, 3)-a(2, 3)*a(3, 2)) * t1
      b (1, 2) = (a(1, 3)*a(3, 2)-a(1, 2)*a(3, 3)) * t1
      b (1, 3) = (a(1, 2)*a(2, 3)-a(1, 3)*a(2, 2)) * t1
      b (2, 1) = (a(2, 3)*a(3, 1)-a(2, 1)*a(3, 3)) * t1
      b (2, 2) = (a(1, 1)*a(3, 3)-a(1, 3)*a(3, 1)) * t1
      b (2, 3) = (a(1, 3)*a(2, 1)-a(1, 1)*a(2, 3)) * t1
      b (3, 1) = (a(2, 1)*a(3, 2)-a(2, 2)*a(3, 1)) * t1
      b (3, 2) = (a(1, 2)*a(3, 1)-a(1, 1)*a(3, 2)) * t1
      b (3, 3) = (a(1, 1)*a(2, 2)-a(1, 2)*a(2, 1)) * t1
      Return
End Subroutine
!EOC
