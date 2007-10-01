
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrfir
! !INTERFACE:
subroutine symrfir(sym,rfir,srfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   sym   : symmetry (in,integer(3,3))
!   rfir  : input interstitial function (in,real(ngrid(1)*ngrid(2)*ngrid(3)))
!   srfir : output interstitial function (out,real(ngrid(1)*ngrid(2)*ngrid(3)))
! !DESCRIPTION:
!   Applies a symmetry to a real intersitial function, $f$, defined on an
!   integer grid. In other words, for each integer vector ${\bf v}$ the routine
!   returns $$ Sf({\bf v})\equiv f(S^{-1}{\bf v}), $$
!   where $S$ is a $3\times 3$ symmetry matrix. Note that the arrays {\tt rfir}
!   and {\tt srfir} represent the functions $f$ and $Sf$ respectively, but
!   stored in the usual way with indices begining at 1.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: sym(3,3)
real(8), intent(in) :: rfir(*)
real(8), intent(out) :: srfir(*)
! local variables
integer i1,i2,i3,iv(3),ir,irs,i,n
real(8) s(3,3),vn(3),t1
real(8) v1(3),v2(3),v3(3)
s(:,:)=dble(sym(:,:))
vn(:)=1.d0/dble(ngrid(:))
ir=0
do i3=0,ngrid(3)-1
  t1=dble(i3)*vn(3)
  v3(:)=t1*s(:,3)
  do i2=0,ngrid(2)-1
    t1=dble(i2)*vn(2)
    v2(:)=t1*s(:,2)+v3(:)
    do i1=0,ngrid(1)-1
      t1=dble(i1)*vn(1)
      v1(:)=t1*s(:,1)+v2(:)
      do i=1,3
        n=ngrid(i)
        iv(i)=modulo(nint(v1(i)*n),n)
      end do
      ir=ir+1
      irs=(iv(3)*ngrid(2)+iv(2))*ngrid(1)+iv(1)+1
      srfir(irs)=rfir(ir)
    end do
  end do
end do
return
end subroutine
!EOC

