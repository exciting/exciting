
! Copyright (C) 2003-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: i3minv
! !INTERFACE:
subroutine i3minv(a,b)
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
implicit none
! arguments
integer, intent(in) :: a(3,3)
integer, intent(out) :: b(3,3)
! local variables
integer m
m=a(1,2)*a(2,3)*a(3,1)-a(1,3)*a(2,2)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
 -a(1,1)*a(2,3)*a(3,2)+a(1,1)*a(2,2)*a(3,3)-a(1,2)*a(2,1)*a(3,3)
if ((m.ne.1).and.(m.ne.-1)) then
  write(*,*)
  write(*,'("Error(i3minv): cannot invert matrix")')
  write(*,'(" Determinant : ",I8)') m
  write(*,*)
  stop
end if
b(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))*m
b(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))*m
b(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))*m
b(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))*m
b(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))*m
b(2,3)=(a(1,3)*a(2,1)-a(1,1)*a(2,3))*m
b(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))*m
b(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))*m
b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))*m
return
end subroutine
!EOC
