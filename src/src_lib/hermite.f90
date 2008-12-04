
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: hermite
! !INTERFACE:
real(8) function hermite(n,x)
! !INPUT/OUTPUT PARAMETERS:
!   n : order of Hermite polynomial (in,integer)
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the $n$th Hermite polynomial. The recurrence relation
!   $$ H_i(x)=2xH_{i-1}(x)-2nH_{i-2}(x), $$
!   with $H_0=1$ and $H_1=2x$, is used. This procedure is numerically stable
!   and accurate to near machine precision for $n\le 20$.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x
! local variables
integer i
real(8) h1,h2,ht
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(hermite): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.gt.20) then
  write(*,*)
  write(*,'("Error(hermite): n out of range : ",I8)') n
  write(*,*)
  stop
end if
if (abs(x).gt.1.d15) then
  write(*,*)
  write(*,'("Error(hermite): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
if (n.eq.0) then
  hermite=1.d0
  return
end if
if (n.eq.1) then
  hermite=2.d0*x
  return
end if
h1=2.d0*x
h2=1.d0
do i=2,n
  ht=2.d0*x*h1-2.d0*dble(i-1)*h2
  h2=h1
  h1=ht
end do
hermite=h1
return
end function
!EOC
