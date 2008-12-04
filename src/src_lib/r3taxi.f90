
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3taxi
! !INTERFACE:
real(8) function r3taxi(x,y)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector 1 (in,real(3))
!   y : input vector 2 (in,real(3))
! !DESCRIPTION:
!   Returns the taxi-cab distance between two real 3-vectors:
!   $d=|x_1-y_1|+|x_2-y_2|+|x_3-y_3|$.
!
! !REVISION HISTORY:
!   Created March 2006 (JKD)
!EOP
!BOC
implicit none
real(8), intent(in) :: x(3)
real(8), intent(in) :: y(3)
r3taxi=abs(x(1)-y(1))+abs(x(2)-y(2))+abs(x(3)-y(3))
return
end function
!EOC

