
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sdelta_fd
! !INTERFACE:
real(8) function sdelta_fd(x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the Fermi-Dirac approximation to the Dirac delta function
!   $$ \tilde\delta(x)=\frac{e^{-x}}{(1+e^{-x})^2}. $$
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x
! local variables
real(8) t1
if (abs(x).gt.50.d0) then
  sdelta_fd=0.d0
  return
end if
t1=exp(-x)
sdelta_fd=t1/((1.d0+t1)**2)
return
end function
!EOC
