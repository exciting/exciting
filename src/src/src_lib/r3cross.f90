!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: r3cross
! !INTERFACE:
!
!
Subroutine r3cross (x, y, z)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector 1 (in,real(3))
!   y : input vector 2 (in,real(3))
!   z : output cross-product (out,real(3))
! !DESCRIPTION:
!   Returns the cross product of two real 3-vectors.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: x (3)
      Real (8), Intent (In) :: y (3)
      Real (8), Intent (Out) :: z (3)
      z (1) = x (2) * y (3) - x (3) * y (2)
      z (2) = x (3) * y (1) - x (1) * y (3)
      z (3) = x (1) * y (2) - x (2) * y (1)
      Return
End Subroutine
!EOC
