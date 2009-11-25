!
!
!
! Copyright (C) 2003-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: i3minv
! !INTERFACE:
!
!
Subroutine i3minv (a, b)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,integer(3,3))
!   b : output matrix (in,integer(3,3))
! !DESCRIPTION:
!   Computes the inverse of a integer $3\times 3$ matrix: $B=A^{-1}$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: a (3, 3)
      Integer, Intent (Out) :: b (3, 3)
! local variables
      Integer :: m
      m = a (1, 2) * a (2, 3) * a (3, 1) - a (1, 3) * a (2, 2) * a (3, &
     & 1) + a (1, 3) * a (2, 1) * a (3, 2) - a (1, 1) * a (2, 3) * a &
     & (3, 2) + a (1, 1) * a (2, 2) * a (3, 3) - a (1, 2) * a (2, 1) * &
     & a (3, 3)
      If ((m .Ne. 1) .And. (m .Ne.-1)) Then
         Write (*,*)
         Write (*, '("Error(i3minv): cannot invert matrix")')
         Write (*, '(" Determinant : ", I8)') m
         Write (*,*)
         Stop
      End If
      b (1, 1) = (a(2, 2)*a(3, 3)-a(2, 3)*a(3, 2)) * m
      b (1, 2) = (a(1, 3)*a(3, 2)-a(1, 2)*a(3, 3)) * m
      b (1, 3) = (a(1, 2)*a(2, 3)-a(1, 3)*a(2, 2)) * m
      b (2, 1) = (a(2, 3)*a(3, 1)-a(2, 1)*a(3, 3)) * m
      b (2, 2) = (a(1, 1)*a(3, 3)-a(1, 3)*a(3, 1)) * m
      b (2, 3) = (a(1, 3)*a(2, 1)-a(1, 1)*a(2, 3)) * m
      b (3, 1) = (a(2, 1)*a(3, 2)-a(2, 2)*a(3, 1)) * m
      b (3, 2) = (a(1, 2)*a(3, 1)-a(1, 1)*a(3, 2)) * m
      b (3, 3) = (a(1, 1)*a(2, 2)-a(1, 2)*a(2, 1)) * m
      Return
End Subroutine
!EOC
