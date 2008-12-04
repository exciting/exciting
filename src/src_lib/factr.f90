
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: factr
! !INTERFACE:
real(8) function factr(n,d)
! !INPUT/OUTPUT PARAMETERS:
!   n : numerator (in,integer)
!   d : denominator (in,integer)
! !DESCRIPTION:
!   Returns the ratio $n!/d!$ for $n,d\ge 0$. Performs no under- or overflow
!   checking.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
integer, intent(in) :: d
! local variables
integer i
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(factr): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
if (d.lt.0) then
  write(*,*)
  write(*,'("Error(factr): d < 0 : ",I8)') d
  write(*,*)
  stop
end if
if (n.lt.d) then
  factr=1.d0/dble(n+1)
  do i=n+2,d
    factr=factr/dble(i)
  end do
else if (n.eq.d) then
  factr=1.d0
else
  factr=dble(d+1)
  do i=d+2,n
    factr=factr*dble(i)
  end do
end if
return
end function
!EOC
