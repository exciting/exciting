!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: r3mtv
! !INTERFACE:
!
!
Subroutine r3mtv (a, x, y)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,real(3,3))
!   x : input vector (in,real(3))
!   y : output vector (out,real(3))
! !DESCRIPTION:
!   Multiplies the transpose of a real $3\times 3$ matrix with a vector.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: a (3, 3)
      Real (8), Intent (In) :: x (3)
      Real (8), Intent (Out) :: y (3)
      y (1) = a (1, 1) * x (1) + a (2, 1) * x (2) + a (3, 1) * x (3)
      y (2) = a (1, 2) * x (1) + a (2, 2) * x (2) + a (3, 2) * x (3)
      y (3) = a (1, 3) * x (1) + a (2, 3) * x (2) + a (3, 3) * x (3)
      Return
End Subroutine
!EOC
