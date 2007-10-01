
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: ztorflm
! !INTERFACE:
subroutine ztorflm(lmax,zflm,rflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   zflm : coefficients of complex spherical harmonic expansion
!          (in,complex((lmax+1)**2)))
!   rflm : coefficients of real spherical harmonic expansion
!          (out,real((lmax+1)**2)))
! !DESCRIPTION:
!   Converts a real function, $z_{lm}$, expanded in terms of complex spherical
!   harmonics into a real spherical harmonic expansion, $r_{lm}$:
!   $$ r_{lm}=\begin{cases}\frac{1}{\sqrt{2}}\Re(z_{lm}+(-1)^m z_{l-m}) & m>0 \\
!    \frac{1}{\sqrt{2}}\Im(-z_{lm}+(-1)^m z_{l-m}) & m<0 \\
!    \Re(z_{lm}) & m=0 \end{cases}\;. $$
!   See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
complex(8), intent(in) :: zflm(*)
real(8), intent(out) :: rflm(*)
! local variables
integer l,m,lm1,lm2
! real constant 1/sqrt(2)
real(8), parameter :: c1=0.7071067811865475244d0
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(ztorflm): lmax < 0 : ",I8)') lmax
  write(*,*)
  stop
end if
do l=0,lmax
  lm1=l**2
  lm2=(l+1)**2+1
  do m=-l,l
    lm1=lm1+1
    lm2=lm2-1
    if (m.gt.0) then
      if (mod(m,2).ne.0) then
        rflm(lm1)=c1*dble(zflm(lm1)-zflm(lm2))
      else
        rflm(lm1)=c1*dble(zflm(lm1)+zflm(lm2))
      end if
    else if (m.lt.0) then
      if (mod(m,2).ne.0) then
        rflm(lm1)=c1*aimag(-zflm(lm1)-zflm(lm2))
      else
        rflm(lm1)=c1*aimag(-zflm(lm1)+zflm(lm2))
      end if
    else
      rflm(lm1)=dble(zflm(lm1))
    end if
  end do
end do
return
end subroutine
!EOC
