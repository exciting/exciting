
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: stheta_fd
! !INTERFACE:
real(8) function stheta_fd(x)
! !INPUT/OUTPUT PARAMETERS:
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the Fermi-Dirac approximation to the Heaviside step function
!   $$ \tilde\Theta(x)=\frac{1}{1+e^{-x}}. $$
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x
if (x.gt.50.d0) then
  stheta_fd=1.d0
  return
end if
if (x.lt.-50.d0) then
  stheta_fd=0.d0
  return
end if
stheta_fd=1.d0/(1.d0+exp(-x))
return
end function
!EOC
