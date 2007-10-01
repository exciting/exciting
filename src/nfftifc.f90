
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: nfftifc
! !INTERFACE:
subroutine nfftifc(n)
! !INPUT/OUTPUT PARAMETERS:
!   n : required/avalable grid size (in,integer)
! !DESCRIPTION:
!   Interface to the grid requirements of the fast Fourier transform routine.
!   Most routines restrict $n$ to specific prime factorisations. This routine
!   returns the next largest grid size allowed by the FFT routine.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(inout) :: n
! local variables
integer i,j
! currently we use primes 2, 3 and 5
integer, parameter :: np=3
integer p(np)
data p / 2,3,5 /
if (n.le.0) then
  write(*,*)
  write(*,'("Error(nfftifc): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
10 continue
i=n
do j=1,np
  do while(mod(i,p(j)).eq.0)
    i=i/p(j)
  end do
end do
if (i.ne.1) then
  n=n+1
  goto 10
end if
end subroutine
!EOC
