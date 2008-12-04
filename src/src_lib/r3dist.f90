
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3dist
! !INTERFACE:
real(8) function r3dist(x,y)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector 1 (in,real(3))
!   y : input vector 2 (in,real(3))
! !DESCRIPTION:
!   Returns the distance between two real 3-vectors: $d=|{\bf x}-{\bf y}|$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x(3)
real(8), intent(in) :: y(3)
r3dist=sqrt((x(1)-y(1))**2+(x(2)-y(2))**2+(x(3)-y(3))**2)
return
end function
!EOC
