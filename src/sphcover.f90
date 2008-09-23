
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sphcover
! !INTERFACE:
subroutine sphcover(n,tp)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of required points (in,integer)
!   tp : (theta, phi) coordinates (out,real(2,n))
! !DESCRIPTION:
!   Produces a set of $N$ points which cover the unit sphere nearly optimally.
!   The points in $(\theta,\phi)$ coordinates are generated using the explicit
!   formula
!   \begin{align*}
!    \theta_k&=\arccos(h_k), \qquad h_k=\frac{2(k-1)}{N-1}-1, \qquad
!    1\le k \le N \\
!    \phi_k&=\left(\phi_{k-1}+C/\sqrt{N(1-h_k^2)}\right)({\rm mod}\;2\pi),
!    \qquad 2\le k\le N-1, \qquad \phi_1=\phi_N=0,
!   \end{align*}
!   where $C=(8\pi/\sqrt{3})^{1/2}$. See E. B. Saff and A. B. J. Kuijlaars,
!   {\it Math. Intell.} {\bf 19}, 5 (1997).
!
! !REVISION HISTORY:
!   Created April 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(out) :: tp(2,n)
! local variables
integer k
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8) c,h,t1,t2
if (n.le.0) then
  write(*,*)
  write(*,'("Error(sphcover): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
c=sqrt(8.d0*pi/sqrt(3.d0))
t1=c/sqrt(dble(n))
tp(1,1)=pi
tp(2,1)=0.d0
do k=2,n-1
  h=dble(2*(k-1))/dble(n-1)-1.d0
  tp(1,k)=acos(h)
  t2=tp(2,k-1)+t1/sqrt(1.d0-h**2)
  tp(2,k)=mod(t2,twopi)
end do
tp(1,n)=0.d0
tp(2,n)=0.d0
return
end subroutine
!EOC
