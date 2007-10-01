
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3dot
! !INTERFACE:
real(8) function r3dot(x,y)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector 1 (in,real(3))
!   y : input vector 2 (in,real(3))
! !DESCRIPTION:
!   Returns the dot-product of two real 3-vectors.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x(3)
real(8), intent(in) :: y(3)
r3dot=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
return
end function
!EOC

