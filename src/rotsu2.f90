
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rotsu2
subroutine rotsu2(rot,det,su2)
! !INPUT/OUTPUT PARAMETERS:
!   rot : input O(3) rotation matrix (in,real(3,3))
!   det : determinant of rotation matrix (out,real)
!   su2 : SU(2) representation of input matrix (out,complex(2,2))
! !DESCRIPTION:
!   Finds the complex ${\rm SU}(2)$ representation of a real ${\rm O}(3)$
!   rotation matrix, which may then be used for rotating 2-spinors. This is done
!   by first determining the rotation matrix in axis-angle coordinates,
!   $(\hat{\bf v},\theta)$, and then using the explicit formula for the
!   spinor rotation matrix:
!   $$ R^{1/2}(\hat{\bf v},\theta)=I\cos\frac{\theta}{2}
!    +i(\hat{\bf v}\cdot\vec{\sigma})\sin\frac{\theta}{2}. $$
!   Note that only the proper part of the rotation matrix is used (i.e.
!   inversion is factored out). See routine {\tt rotaxang}.
!
! !REVISION HISTORY:
!   Created August 2007 (JKD)
!   Reversed rotation direction, February 2008 (L. Nordstrom)
!EOP
!BOC
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
su2(1,1)=cmplx(cs,v(3)*sn,8)
su2(1,2)=cmplx(v(2)*sn,v(1)*sn,8)
su2(2,1)=cmplx(-v(2)*sn,v(1)*sn,8)
su2(2,2)=cmplx(cs,-v(3)*sn,8)
return
end subroutine
!EOC

