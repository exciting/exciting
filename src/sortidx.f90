
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sortidx
! !INTERFACE:
subroutine sortidx(n,a,idx)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of elements in array (in,integer)
!   idx : permutation index (out,integer(n))
!   a   : real array (in,real(n))
! !DESCRIPTION:
!   Finds the permutation index {\tt idx} which sorts the real array {\tt a}
!   into ascending order. No sorting of the array {\tt a} itself is performed.
!   Uses the heapsort algorthim.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Included tolerance eps, April 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: a(n)
integer, intent(out) :: idx(n)
! local variables
integer i,j,k,l,m
! tolerance for deciding if one number is smaller than another
real(8), parameter :: eps=1.d-14
if (n.le.0) then
  write(*,*)
  write(*,'("Error(sortidx): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
do i=1,n
  idx(i)=i
end do
if (n.eq.1) return
l=n/2+1
k=n
10 continue
if (l.gt.1) then
  l=l-1
  m=idx(l)
else
  m=idx(k)
  idx(k)=idx(1)
  k=k-1
  if (k.eq.1) then
    idx(1)=m
    return
  end if
end if
i=l
j=l+l
20 continue
if (j.le.k) then
  if (j.lt.k) then
    if (a(idx(j)).lt.a(idx(j+1))+eps) j=j+1
  end if
  if (a(m).lt.a(idx(j))+eps) then
    idx(i)=idx(j)
    i=j
    j=j+j
  else
    j=k+1
  end if
  goto 20
end if
idx(i)=m
goto 10
end subroutine
!EOC
