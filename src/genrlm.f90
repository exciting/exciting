
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genrlm
! !INTERFACE:
subroutine genrlm(lmax,tp,rlm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   tp   : (theta, phi) coordinates (in,real(2))
!   rlm  : array of real spherical harmonics (out,real((lmax+1)**2))
! !DESCRIPTION:
!   Generates a sequence of real spherical harmonics evaluated at angles
!   $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values are returned in a packed
!   array {\tt rlm} indexed with $j=l(l+1)+m+1$. Real spherical harmonics are
!   defined by
!   $$ R_{lm}(\theta,\phi)= \begin{cases}
!     \sqrt{2}\,\Re\{Y_{lm}(\theta,\phi)\} & m>0 \\
!     \sqrt{2}\,\Im\{Y_{lm}(\theta,\phi)\} & m<0 \\
!     \Re\{Y_{lm}(\theta,\phi)\} & m=0
!    \end{cases}, $$
!   where $Y_{lm}$ are the complex spherical harmonics. These functions are
!   orthonormal and complete and may be used for expanding real-valued functions
!   on the sphere. This routine is numerically stable and accurate to near
!   machine precision for $l\le 50$. See routine {\tt genylm}.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: tp(2)
real(8), intent(out) :: rlm(*)
! local variables
integer lmmax,l,m,lm
real(8), parameter :: sqtwo=1.4142135623730950488d0
! allocatable arrays
complex(8), allocatable :: ylm(:)
if ((lmax.lt.0).or.(lmax.gt.50)) then
  write(*,*)
  write(*,'("Error(genrlm): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
allocate(ylm(lmmax))
! generate complex spherical harmonics
call genylm(lmax,tp,ylm)
! convert to real spherical harmonics
lm=0
do l=0,lmax
  do m=-l,-1
    lm=lm+1
    rlm(lm)=sqtwo*aimag(ylm(lm))
  end do
  lm=lm+1
  rlm(lm)=dble(ylm(lm))
  do m=1,l
    lm=lm+1
    rlm(lm)=sqtwo*dble(ylm(lm))
  end do
end do
deallocate(ylm)
return
end subroutine
!EOC
