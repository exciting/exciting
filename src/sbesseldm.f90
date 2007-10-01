
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sbesseldm
! !INTERFACE:
subroutine sbesseldm(m,lmax,x,djl)
! !INPUT/OUTPUT PARAMETERS:
!   m    : order of derivatve (in,integer)
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   djl  : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes the $m$th derivative of the spherical Bessel function of the first
!   kind, $j_l(x)$, for argument $x$ and $l=0,1,\ldots,l_{\rm max}$. For
!   $x\ge 1$ this is done by repeatedly using the relations
!   \begin{align*}
!    \frac{d}{dx}j_l(x)&=\frac{l}{x}j_l(x)-j_{l+1}(x) \\
!    j_{l+1}(x)&=\frac{2l+1}{x}j_l(x)-j_{l-1}(x).
!   \end{align*}
!   While for $x<1$ the series expansion of the Bessel function is used
!   $$ \frac{d^m}{dx^m}j_l(x)=\sum_{i=0}^{\infty}
!    \frac{(2i+l)!}{(-2)^ii!(2i+l-m)!(2i+2l+1)!!}x^{2i+l-m}. $$
!   This procedure is numerically stable and accurate to near machine precision
!   for $l\le 30$ and $m\le 6$.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Modified to return an array of values, October 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m
integer, intent(in) :: lmax
real(8), intent(in) :: x
real(8), intent(out) :: djl(0:lmax)
! local variables
integer, parameter :: maxm=6
integer, parameter :: maxns=20
integer i,j,l,i0
real(8) t1,sum,x2
integer a(0:maxm+1),a1(0:maxm+1)
integer b(0:maxm+1),b1(0:maxm+1)
! automatic arrays
real(8) jl(0:lmax+1)
! external functions
real(8) factnm,factr
external factnm,factr
if ((m.lt.0).or.(m.gt.maxm)) then
  write(*,*)
  write(*,'("Error(sbesseldm): m out of range : ",I8)') m
  write(*,*)
  stop
end if
if ((lmax.lt.0).or.(lmax.gt.30)) then
  write(*,*)
  write(*,'("Error(sbesseldm): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
if ((x.lt.0.d0).or.(x.gt.1.d5)) then
  write(*,*)
  write(*,'("Error(sbesseldm): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
if (m.eq.0) then
  call sbessel(lmax,x,djl)
  return
end if
if (x.gt.1.d0) then
  call sbessel(lmax+1,x,jl)
  do l=0,lmax
    a(1:m+1)=0
    a(0)=1
    a1(0:m+1)=0
    do i=1,m
      b(0)=0
      b1(0)=0
      do j=0,i
        b(j+1)=a(j)*(l-j)
        b1(j+1)=-a1(j)*(j+l+2)
      end do
      do j=0,i
        b1(j)=b1(j)-a(j)
        b(j)=b(j)+a1(j)
      end do
      a(0:i+1)=b(0:i+1)
      a1(0:i+1)=b1(0:i+1)
    end do
    t1=1.d0
    sum=dble(a(0))*jl(l)+dble(a1(0))*jl(l+1)
    do i=1,m+1
      t1=t1*x
      sum=sum+(dble(a(i))*jl(l)+dble(a1(i))*jl(l+1))/t1
    end do
    djl(l)=sum
  end do
else
  x2=x**2
  do l=0,lmax
    i0=max((m-l+1)/2,0)
    j=2*i0+l-m
    if (j.eq.0) then
      t1=1.d0
    else
      t1=x**j
    end if
    t1=factr(j+m,j)*t1/(factnm(i0,1)*factnm(j+l+m+1,2)*dble((-2)**i0))
    sum=t1
    do i=i0+1,maxns
      j=2*i+l
      t1=-t1*dble((j-1)*j)*x2/dble((j-l)*(j-m-1)*(j-m)*(j+l+1))
      if (abs(t1).le.1.d-40) goto 10
      sum=sum+t1
    end do
10 continue
    djl(l)=sum
  end do
end if
return
end subroutine
!EOC
