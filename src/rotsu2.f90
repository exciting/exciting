
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rotsu2(rot,det,su2)
implicit none
! arguments
real(8), intent(in) :: rot(3,3)
real(8), intent(out) :: det
complex(8), intent(out) :: su2(2,2)
! local variables
real(8), parameter :: eps=1.d-6
real(8) v(3),th,cs,sn
! determine the axis and angle of rotation for the rotation matrix
call rotaxang(eps,rot,det,v,th)
! determine the SU(2) representation of the rotation matrix
cs=cos(th/2.d0)
sn=sin(th/2.d0)
su2(1,1)=cmplx(cs,-v(3)*sn,8)
su2(1,2)=cmplx(-v(2)*sn,-v(1)*sn,8)
su2(2,1)=cmplx(v(2)*sn,-v(1)*sn,8)
su2(2,2)=cmplx(cs,v(3)*sn,8)
return
end subroutine

