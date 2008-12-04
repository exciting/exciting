
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3mdet
! !INTERFACE:
real(8) function r3mdet(a)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,real(3,3))
! !DESCRIPTION:
!   Returns the determinant of a real $3\times 3$ matrix $A$.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3)
r3mdet=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
      +a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
      +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
return
end function
!EOC
