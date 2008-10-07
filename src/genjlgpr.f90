
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjlgpr(lmax,gpc,jlgpr)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: gpc(ngvec)
real(8), intent(out) :: jlgpr(0:lmax,ngvec,nspecies)
! local variables
integer is,ig
real(8) t1
do is=1,nspecies
  do ig=1,ngvec
    t1=gpc(ig)*rmt(is)
    call sbessel(lmax,t1,jlgpr(:,ig,is))
  end do
end do
return
end subroutine

