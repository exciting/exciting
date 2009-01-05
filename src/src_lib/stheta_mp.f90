
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: stheta_mp
! !INTERFACE:
real(8) function stheta_mp(n,x)
! !INPUT/OUTPUT PARAMETERS:
!   n : order (in,integer)
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the smooth approximation to the Heaviside step function of order
!   $N$ given by Methfessel and Paxton, {\it Phys. Rev. B} {\bf 40}, 3616
!   (1989),
!   $$ \tilde\Theta(x)=1-S_N(x) $$
!   where
!   \begin{align*}
!    S_N(x)&=S_0(x)+\sum_{i=1}^N \frac{(-1)^i}{i!4^n\sqrt\pi} H_{2i-1}(x)
!     e^{-x^2},\\
!    S_0(x)&=\frac{1}{2}(1-{\rm erf}(x))
!   \end{align*}
!   and $H_j$ is the $j$th-order Hermite polynomial. This procedure is numerically
!   stable and accurate to near machine precision for $N\le 10$.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x
! local variables
integer i
real(8), parameter :: sqpi=1.7724538509055160273d0
real(8) sum,t1
! external functions
real(8) factnm,hermite,erf
external factnm,hermite,erf
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(stheta_mp): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.gt.10) then
  write(*,*)
  write(*,'("Error(stheta_mp): n out of range : ",I8)') n
  write(*,*)
  stop
end if
if (x.lt.-12.d0) then
  stheta_mp=0.d0
  return
end if
if (x.gt.12.d0) then
  stheta_mp=1.d0
  return
end if
sum=0.5d0*(1.d0-erf(x))
do i=1,n
  t1=1.d0/(factnm(i,1)*dble(4**i)*sqpi)
  if (mod(i,2).ne.0) t1=-t1
  sum=sum+t1*hermite(2*i-1,x)*exp(-x**2)
end do
stheta_mp=1.d0-sum
return
end function
!EOC
