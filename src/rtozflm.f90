
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rtozflm
! !INTERFACE:
subroutine rtozflm(lmax,rflm,zflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   rflm : coefficients of real spherical harmonic expansion
!          (in,real((lmax+1)**2)))
!   zflm : coefficients of complex spherical harmonic expansion
!          (out,complex((lmax+1)**2)))
! !DESCRIPTION:
!   Converts a real function, $r_{lm}$, expanded in terms of real spherical
!   harmonics into a complex spherical harmonic expansion, $z_{lm}$:
!   $$ z_{lm}=\begin{cases} \frac{1}{\sqrt{2}}(r_{lm}+i(-1)^mr_{l-m}) & m>0 \\
!    \frac{1}{\sqrt{2}}((-1)^mr_{l-m}-ir_{lm}) & m<0 \\
!    r_{lm} & m=0 \end{cases}\;. $$
!   See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: rflm(*)
complex(8), intent(out) :: zflm(*)
! local variables
integer l,m,lm1,lm2
! real constant 1/sqrt(2)
real(8), parameter :: c1=0.7071067811865475244d0
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(rtozflm): lmax < 0 : ",I8)') lmax
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
        zflm(lm1)=c1*cmplx(rflm(lm1),-rflm(lm2),8)
      else
        zflm(lm1)=c1*cmplx(rflm(lm1),rflm(lm2),8)
      end if
    else if (m.lt.0) then
      if (mod(m,2).ne.0) then
        zflm(lm1)=c1*cmplx(-rflm(lm2),-rflm(lm1),8)
      else
        zflm(lm1)=c1*cmplx(rflm(lm2),-rflm(lm1),8)
      end if
    else
      zflm(lm1)=rflm(lm1)
    end if
  end do
end do
return
end subroutine
!EOC
