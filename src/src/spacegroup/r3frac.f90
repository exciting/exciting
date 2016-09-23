
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3frac
! !INTERFACE:
subroutine r3frac(eps,v,iv)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
!   iv  : integer parts of v (out,integer(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!   The integer components of {\tt v} are returned in the variable {\tt iv}.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: v(3)
integer, intent(out) :: iv(3)
! local variables
integer i
do i=1,3
  iv(i)=int(v(i))
  v(i)=v(i)-dble(iv(i))
  if (v(i).lt.0.d0) then
    v(i)=v(i)+1.d0
    iv(i)=iv(i)-1
  end if
  if (1.d0-v(i).lt.eps) then
    v(i)=0.d0
    iv(i)=iv(i)+1
  end if
  if (v(i).lt.eps) then
    v(i)=0.d0
  end if
end do
return
end subroutine
!EOC
