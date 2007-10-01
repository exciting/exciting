
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findigp0(ngp,gpc,igp0)
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: gpc(ngp)
integer, intent(out) :: igp0
! local variables
integer igp
real(8), parameter :: eps=1.d-14
real(8) t1
igp0=1
t1=gpc(igp0)+eps
do igp=2,ngp
  if (gpc(igp).lt.t1) then
    igp0=igp
    t1=gpc(igp)+eps
  end if
end do
return
end subroutine

