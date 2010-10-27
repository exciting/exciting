
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: spline
! !INTERFACE:
subroutine spline4(n,x,ld,f,cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   ld : leading dimension (in,integer)
!   f  : input data array (in,real(ld,n))
!   cf : cubic spline coefficients (1,2,3) and work space (4) (out,real(4,n))
! !DESCRIPTION:
!   Calculates the coefficients of a cubic spline fitted to input data. In other
!   words, given a set of data points $f_i$ defined at $x_i$, where
!   $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!   $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!   is the interpolating function for $x\in[x_i,x_{i+1})$. This is done by
!   determining the end-point coefficients $c_2^1$ and $c_2^n$ from the first
!   and last three points, and then solving the tridiagonal system
!   $$ d_{i-1}c_2^{i-1}+2(d_{i-1}+d_i)c_2^i+d_ic_2^{i+1}
!    =3\left(\frac{f_{i+1}-f_i}{d_i}-\frac{f_i-f_{i-1}}{d_{i-1}}\right), $$
!   where $d_i=x_{i+1}-x_i$, for the intermediate coefficients.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Improved speed and accuracy, April 2006 (JKD)
!   Optimisations and improved end-point coefficients, February 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
integer, intent(in) :: ld
real(8), intent(in) :: f(ld,n)
real(8), intent(out) :: cf(4,n)
! local variables
integer i
real(8) t1,t2,t3,t4
if (n.le.0) then
  write(*,*)
  write(*,'("Error(spline): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.eq.1) then
  cf(:,1)=0.d0
  return
end if
if (n.eq.2) then
  cf(1,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
  cf(2:3,1)=0.d0
  cf(1,2)=cf(1,1)
  cf(2:3,2)=0.d0
  return
end if
cf(4,1)=1.d0/(x(2)-x(1))
cf(1,1)=cf(4,1)*(f(1,2)-f(1,1))
cf(4,2)=1.d0/(x(3)-x(2))
cf(1,2)=cf(4,2)*(f(1,3)-f(1,2))
cf(2,1)=1.d0
! estimate second derivative at the first point
cf(3,1)=(cf(1,2)-cf(1,1))/(x(3)-x(1))
! use Gaussian elimination to solve tridiagonal system
t1=(x(2)-x(1))*cf(4,2)
t2=t1*cf(2,1)
t3=1.d0/(2.d0*(t1+1.d0))
cf(2,2)=t3
t4=3.d0*(cf(1,2)-cf(1,1))*cf(4,2)-t2*cf(3,1)
cf(3,2)=t4
do i=3,n-1
  cf(4,i)=1.d0/(x(i+1)-x(i))
  cf(1,i)=cf(4,i)*(f(1,i+1)-f(1,i))
  t1=(x(i)-x(i-1))*cf(4,i)
  t2=t1*t3
  t3=1.d0/(2.d0*t1+2.d0-t2)
  cf(2,i)=t3
  t4=3.d0*(cf(1,i)-cf(1,i-1))*cf(4,i)-t2*t4
  cf(3,i)=t4
end do
! estimate second derivative at the last point
t3=(cf(1,n-1)-cf(1,n-2))/(x(n)-x(n-2))
cf(3,n)=t3
cf(2,n)=t3
t3=cf(3,n-1)-t3
cf(3,n-1)=t3
t2=cf(2,n-1)
cf(2,n-1)=t2*t3
do i=n-2,2,-1
  t3=cf(3,i)-t2*t3
  cf(3,i)=t3
  t2=cf(2,i)
  cf(2,i)=t2*t3
end do
cf(2,1)=cf(2,1)*cf(3,1)
do i=1,n-1
  t1=0.3333333333333333333d0*(cf(2,i+1)-cf(2,i))
  cf(3,i)=t1*cf(4,i)
  cf(1,i)=cf(1,i)-(cf(2,i)+t1)*(x(i+1)-x(i))
end do
! determine end-point coefficients
t1=x(n)-x(n-1)
cf(1,n)=cf(1,n-1)+(2.d0*cf(2,n-1)+3.d0*cf(3,n-1)*t1)*t1
cf(3,n)=cf(3,n-1)
return
end subroutine
!EOC

