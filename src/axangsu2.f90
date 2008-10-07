
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: axangsu2
subroutine axangsu2(v,th,su2)
! !INPUT/OUTPUT PARAMETERS:
!   v   : rotation axis vector (in,real(3))
!   th  : rotation angle (in,real)
!   su2 : SU(2) representation of rotation (out,complex(2,2))
! !DESCRIPTION:
!   Finds the complex ${\rm SU}(2)$ representation of a rotation defined by an
!   axis vector $\hat{\bf v}$ and angle $\theta$. The spinor rotation matrix is
!   given explicitly by
!   $$ R^{1/2}(\hat{\bf v},\theta)=I\cos\frac{\theta}{2}
!    +i(\hat{\bf v}\cdot\vec{\sigma})\sin\frac{\theta}{2}. $$
!
! !REVISION HISTORY:
!   Created August 2007 (JKD)
!   Reversed rotation direction, February 2008 (L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: v(3)
real(8), intent(in) :: th
complex(8), intent(out) :: su2(2,2)
! local variables
real(8), parameter :: eps=1.d-6
real(8) vn(3),cs,sn,t1
t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
if (t1.lt.eps) then
  write(*,*)
  write(*,'("Error(axangsu2): zero length axis vector")')
  write(*,*)
  stop
end if
! normalise the vector
vn(:)=v(:)/t1
cs=cos(th/2.d0)
sn=sin(th/2.d0)
su2(1,1)=cmplx(cs,vn(3)*sn,8)
su2(1,2)=cmplx(vn(2)*sn,vn(1)*sn,8)
su2(2,1)=cmplx(-vn(2)*sn,vn(1)*sn,8)
su2(2,2)=cmplx(cs,-vn(3)*sn,8)
return
end subroutine
!EOC

