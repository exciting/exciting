
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!*************************************************************************
! PURPOSE: Evaluate inverse matrix of size (3x3)  
! DATA:    
! COMMENT: 
!*************************************************************************

subroutine r3minv(a,b)

! arguments
real(8):: a(3,3)
real(8):: b(3,3)
!
! local variables
real(8) c(3,3),t1
t1=a(1,2)*a(2,3)*a(3,1)-a(1,3)*a(2,2)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
  -a(1,1)*a(2,3)*a(3,2)+a(1,1)*a(2,2)*a(3,3)-a(1,2)*a(2,1)*a(3,3)
if (abs(t1).lt.1.d-40) then

  write(*,*)
  write(*,'("Error(r3minv): singular matrix")')
  write(*,*)
  stop
end if
t1=1.d0/t1
c(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))*t1
c(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))*t1
c(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))*t1

c(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))*t1
c(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))*t1
c(2,3)=(a(1,3)*a(2,1)-a(1,1)*a(2,3))*t1
c(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))*t1
c(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))*t1
c(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))*t1

! copy to output matrix
b(:,:)=c(:,:)
return
end subroutine
!*************************************************************************
