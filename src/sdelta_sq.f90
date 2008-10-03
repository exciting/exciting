
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sdelta_sq
! !INTERFACE:
real(8) function sdelta_sq(x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the square-wave pulse approximation to the Dirac delta function
!   $$ \tilde\delta(x)=\left\{\begin{array}{ll}
!    1 & \quad |x|\le 1/2 \\
!    0 & \quad |x|>1/2 \end{array}\right. $$
!
! !REVISION HISTORY:
!   Created July 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x
if (abs(x).le.0.5d0) then
  sdelta_sq=1.d0
else
  sdelta_sq=0.d0
end if
return
end function
!EOC
