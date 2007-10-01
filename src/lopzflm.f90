
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: lopzflm
! !INTERFACE:
subroutine lopzflm(lmax,zflm,ld,zlflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   zflm  : coefficients of input spherical harmonic expansion
!           (in,complex((lmax+1)**2))
!   ld    : leading dimension (in,integer)
!   zlflm : coefficients of output spherical harmonic expansion
!           (out,complex(ld,3))
! !DESCRIPTION:
!   Applies the angular momentum operator ${\bf L}$ to a function expanded in
!   terms of complex spherical harmonics. This makes use of the identities
!   \begin{align*}
!    (L_x+iL_y)Y_{lm}(\theta,\phi)&=\sqrt{(l-m)(l+m+1)}Y_{lm+1}(\theta,\phi)\\
!    (L_x-iL_y)Y_{lm}(\theta,\phi)&=\sqrt{(l+m)(l-m+1)}Y_{lm-1}(\theta,\phi)\\
!    L_zY_{lm}(\theta,\phi)&=mY_{lm}(\theta,\phi).
!   \end{align*}
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
complex(8), intent(in) :: zflm(*)
integer, intent(in) :: ld
complex(8), intent(out) :: zlflm(ld,3)
! local variables
integer l,m,lm
real(8) t1
complex(8) zt1
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(lopzflm): lmax < 0 : ",I8)') lmax
  write(*,*)
  stop
end if
lm=0
do l=0,lmax
  do m=-l,l
    lm=lm+1
    if (m.eq.-l) then
      zlflm(lm,1)=0.d0
      zlflm(lm,2)=0.d0
    end if
    if (m.lt.l) then
      t1=0.5d0*sqrt(dble((l-m)*(l+m+1)))
      zt1=t1*zflm(lm)
      zlflm(lm+1,1)=zt1
      zlflm(lm+1,2)=cmplx(aimag(zt1),-dble(zt1),8)
    end if
    if (m.gt.-l) then
      t1=0.5d0*sqrt(dble((l+m)*(l-m+1)))
      zt1=t1*zflm(lm)
      zlflm(lm-1,1)=zlflm(lm-1,1)+zt1
      zlflm(lm-1,2)=zlflm(lm-1,2)+cmplx(-aimag(zt1),dble(zt1),8)
    end if
    if (m.ne.0) then
      zlflm(lm,3)=dble(m)*zflm(lm)
    else
      zlflm(lm,3)=0.d0
    end if
  end do
end do
return
end subroutine
!EOC
