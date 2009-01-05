
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fderiv
! !INTERFACE:
subroutine fderiv(m,n,x,f,g,cf)
! !INPUT/OUTPUT PARAMETERS:
!   m  : order of derivative (in,integer)
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   f  : function array (in,real(n))
!   g  : (anti-)derivative of f (out,real(n))
!   cf : spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!   routine computes the $m$th derivative of $f$ at each point. If $m<0$ the
!   anti-derivative of $f$ given by
!   $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!   is calculated by fitting the function to a clamped cubic spline. See routine
!   {\tt spline}.
!
! !REVISION HISTORY:
!   Created May 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(in) :: f(n)
real(8), intent(out) :: g(n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) dx
if (n.le.0) then
  write(*,*)
  write(*,'("Error(fderiv): invalid number of points : ",I8)') n
  write(*,*)
  stop
end if
if (m.eq.0) then
  g(:)=f(:)
  return
end if
if (m.ge.4) then
  g(:)=0.d0
  return
end if
! high accuracy (anti-)derivatives from a clamped spline fit to the data
call spline(n,x,1,f,cf)
select case(m)
case(:-1)
  g(1)=0.d0
  do i=1,n-1
    dx=x(i+1)-x(i)
    g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx+0.3333333333333333333d0*cf(2,i))*dx &
     +0.5d0*cf(1,i))*dx+f(i))*dx
  end do
case(1)
  g(:)=cf(1,:)
case(2)
  g(:)=2.d0*cf(2,:)
case(3)
  g(:)=6.d0*cf(3,:)
end select
return
end subroutine
!EOC

