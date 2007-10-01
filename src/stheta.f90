
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: stheta
! !INTERFACE:
real(8) function stheta(stype,x)
! !INPUT/OUTPUT PARAMETERS:
!   stype : smearing type (in,integer)
!   x     : real argument (in,real)
! !DESCRIPTION:
!   Returns the Heaviside step function corresponding to the smooth approximation
!   to the Dirac delta function:
!   $$ \tilde\Theta(x)=\int_{-\infty}^x dt\,\tilde\delta(t). $$
!   See function {\tt sdelta} for details.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: stype
real(8), intent(in) :: x
! external functions
real(8) stheta_mp,stheta_fd
external stheta_mp,stheta_fd
stheta=0.d0
select case(stype)
case(0)
  stheta=stheta_mp(0,x)
  return
case(1)
  stheta=stheta_mp(1,x)
  return
case(2)
  stheta=stheta_mp(2,x)
  return
case(3)
  stheta=stheta_fd(x)
  return
case default
  write(*,*)
  write(*,'("Error(stheta): sytpe not defined : ",I8)') stype
  write(*,*)
  stop
end select
end function
!EOC
