
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: i3mtv
! !INTERFACE:
subroutine i3mtv(a,x,y)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,integer(3,3))
!   x : input vector (in,integer(3))
!   y : output vector (out,integer(3))
! !DESCRIPTION:
!   Multiplies the transpose of an integer $3\times 3$ matrix with a vector.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: a(3,3)
integer, intent(in) :: x(3)
integer, intent(out) :: y(3)
y(1)=a(1,1)*x(1)+a(2,1)*x(2)+a(3,1)*x(3)
y(2)=a(1,2)*x(1)+a(2,2)*x(2)+a(3,2)*x(3)
y(3)=a(1,3)*x(1)+a(2,3)*x(2)+a(3,3)*x(3)
return
end subroutine
!EOC

