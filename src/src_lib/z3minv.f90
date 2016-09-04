!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: z3minv
! !INTERFACE:
!
!
Subroutine z3minv (a, b)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,complex(3,3))
!   b : output matrix (out,complex(3,3))
! !DESCRIPTION:
!   Computes the inverse of a complex $3\times 3$ matrix.
!
! !REVISION HISTORY:
!   Created March 2014 (STK)
!EOP
!BOC
      Implicit None
! arguments
      Complex (8), Intent (In) :: a (3, 3)
      Complex (8), Intent (Out) :: b (3, 3)
! local variables
      Complex (8) :: t1
      t1 = a (1, 2) * a (2, 3) * a (3, 1) - a (1, 3) * a (2, 2) * a (3, &
     & 1) + a (1, 3) * a (2, 1) * a (3, 2) - a (1, 1) * a (2, 3) * a &
     & (3, 2) + a (1, 1) * a (2, 2) * a (3, 3) - a (1, 2) * a (2, 1) * &
     & a (3, 3)
      If (Abs(t1) .Lt. 1.d-40) Then
         Write (*,*)
         Write (*, '("Error(z3minv): singular matrix")')
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
