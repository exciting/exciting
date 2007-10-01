
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gcd
! !INTERFACE:
integer function gcd(x,y)
! !INPUT/OUTPUT PARAMETERS:
!   x : first integer (in,integer)
!   y : second integer (in,integer)
! !DESCRIPTION:
!   Computes the greatest common divisor (GCD) of two integers using Euclid's
!   algorithm.
!
! !REVISION HISTORY:
!   Created September 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: x
integer, intent(in) :: y
! local variables
integer a,b,c
if (x.le.0) then
  write(*,*)
  write(*,'("Error(gcd): x <= 0 : ",I8)') x
  write(*,*)
  stop
end if
if (y.le.0) then
  write(*,*)
  write(*,'("Error(gcd): y <= 0 : ",I8)') y
  write(*,*)
  stop
end if
if (x.ge.y) then
  a=x
  b=y
else
  a=y
  b=x
end if
10 continue
c=mod(a,b)
a=b
b=c
if (c.gt.0) goto 10
gcd=a
return
end function
!EOC
