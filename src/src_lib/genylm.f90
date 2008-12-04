
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genylm
! !INTERFACE:
subroutine genylm(lmax,tp,ylm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   tp   : (theta, phi) coordinates (in,real(2))
!   ylm  : array of spherical harmonics (out,complex((lmax+1)**2))
! !DESCRIPTION:
!   Generates a sequence of spherical harmonics, including the Condon-Shortley
!   phase, evaluated at angles $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values
!   are returned in a packed array {\tt ylm} indexed with $j=l(l+1)+m+1$. The
!   algorithm of Masters and Richards-Dinger is used, {\it Geophys. J. Int.}
!   {\bf 135}, 307 (1998). This routine is numerically stable and accurate to
!   near machine precision for $l\le 50$.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!   Improved stability, December 2005 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: tp(2)
complex(8), intent(out) :: ylm(*)
! local variables
integer l,m,lm1,lm2
real(8), parameter :: fourpi=12.566370614359172954d0
real(8) sn,cs,dx,sum,t1
! automatic arrays
real(8) x(0:lmax)
complex(8) z(lmax)
if ((lmax.lt.0).or.(lmax.gt.50)) then
  write(*,*)
  write(*,'("Error(genylm): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
ylm(1)=0.28209479177387814347d0
if (lmax.eq.0) return
sn=sin(tp(1))
cs=cos(tp(1))
! phase factors exp(i*m*phi)
do m=1,lmax
  t1=dble(m)*tp(2)
  z(m)=cmplx(cos(t1),sin(t1),8)
end do
do l=1,lmax
  if (mod(l,2).eq.0) then
    x(l)=1.d0
  else
    x(l)=-1.d0
  end if
! recursion loop
  dx=0.d0
  do m=l,1,-1
    t1=sqrt(dble((l+m)*(l-m+1)))
    x(m-1)=-(sn*dx+dble(2*m)*cs*x(m))/t1
    dx=sn*x(m)*t1
  end do
! rescale values and multiply with phase factors
  t1=sn
  sum=0.d0
  do m=1,l
    x(m)=t1*x(m)
    sum=sum+x(m)**2
    t1=t1*sn
  end do
  sum=2.d0*sum+x(0)**2
  t1=sqrt(dble(2*l+1)/(fourpi*sum))
  lm1=l*(l+1)+1
  lm2=lm1
  ylm(lm1)=t1*x(0)
  do m=1,l
    lm1=lm1+1
    lm2=lm2-1
    ylm(lm1)=t1*x(m)*z(m)
    ylm(lm2)=conjg(ylm(lm1))
    if (mod(m,2).ne.0) ylm(lm2)=-ylm(lm2)
  end do
end do
return
end subroutine
!EOC
