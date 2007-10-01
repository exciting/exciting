
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: spline
! !INTERFACE:
subroutine spline(n,x,ld,f,cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   ld : leading dimension (in,integer)
!   f  : input data array (in,real(n))
!   cf : cubic spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Calculates the coefficients of a clamped cubic spline fitted to input data.
!   In other words, given a set of data points $f_i$ defined at $x_i$, where
!   $i=1\ldots n$, the coefficients $c_{ij}$ are determined such that
!   $$ y_i(x)=f_i+c_{i1}(x-x_i)+c_{i2}(x-x_i)^2+c_{i3}(x-x_i)^3, $$
!   is the interpolating function for $x\in[x_i,x_{i+1})$.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Improved speed and accuracy, April 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
integer, intent(in) :: ld
real(8), intent(in) :: f(ld,n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) t1,t2
! automatic arrays
real(8) w(n)
do i=1,n-1
  cf(3,i)=1.d0/(x(i+1)-x(i))
  cf(1,i)=cf(3,i)*(f(1,i+1)-f(1,i))
end do
t1=0.5d0/(x(3)-x(1))
cf(2,1)=t1*(cf(1,2)-cf(1,1))
do i=2,n-1
  t2=x(i)-x(i-1)
  w(i)=t2*t1
  t1=1.d0/(2.d0*(x(i+1)-x(i-1))-t2*w(i))
  cf(2,i)=t1*(cf(1,i)-cf(1,i-1)-t2*cf(2,i-1))
end do
! back-substitution
do i=n-2,1,-1
  cf(2,i)=cf(2,i)-w(i+1)*cf(2,i+1)
end do
! determine coefficients
cf(2,n)=0.d0
do i=1,n-1
  t2=cf(2,i+1)-cf(2,i)
  cf(3,i)=t2*cf(3,i)
  cf(2,i)=3.d0*cf(2,i)
  cf(1,i)=cf(1,i)-(cf(2,i)+t2)*(x(i+1)-x(i))
end do
! end-point coefficients for extrapolation
t1=x(n)-x(n-1)
cf(1,n)=(3.d0*cf(3,n-1)*t1+2.d0*cf(2,n-1))*t1+cf(1,n-1)
cf(2,n)=6.d0*cf(3,n-1)*t1+2.d0*cf(2,n-1)
cf(3,n)=cf(3,n-1)
return
end subroutine
!EOC

