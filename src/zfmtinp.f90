
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: zfmtinp
! !INTERFACE:
complex(8) function zfmtinp(lmax,nr,r,ld,zfmt1,zfmt2)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   ld    : leading dimension (in,integer)
!   zfmt1 : first complex function inside muffin-tin (in,complex(ld,nr))
!   zfmt2 : second complex function inside muffin-tin (in,complex(ld,nr))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions in the muffin-tin. In
!   other words, given two complex functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)Y_{lm}
!    (\hat{\bf r}), $$
!   the function returns
!   $$ I=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}\int f_{lm}^{1*}(r)
!    f_{lm}^2(r)r^2\,dr\;. $$
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
integer, intent(in) :: ld
complex(8), intent(in) :: zfmt1(ld,nr)
complex(8), intent(in) :: zfmt2(ld,nr)
! local variables
integer lmmax,ir
real(8) t1,t2
complex(8) zt1
! automatic arrays
real(8) fr1(nr),fr2(nr),gr(nr),cf(3,nr)
! external functions
complex(8) zdotc
external zdotc
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(zfmtinp): lmax < 0 : ",I8)') lmax
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
do ir=1,nr
  t1=r(ir)**2
  zt1=zdotc(lmmax,zfmt1(1,ir),1,zfmt2(1,ir),1)
  fr1(ir)=t1*dble(zt1)
  fr2(ir)=t1*aimag(zt1)
end do
call fderiv(-1,nr,r,fr1,gr,cf)
t1=gr(nr)
call fderiv(-1,nr,r,fr2,gr,cf)
t2=gr(nr)
zfmtinp=cmplx(t1,t2,8)
return
end function
!EOC

