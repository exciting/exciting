
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: reciplat
! !INTERFACE:
subroutine reciplat
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the reciprocal lattice vectors from the real-space lattice vectors
!   \begin{align*}
!     {\bf b}_1&=\frac{2\pi}{s}({\bf a}_2\times{\bf a}_3)\\
!     {\bf b}_2&=\frac{2\pi}{s}({\bf a}_3\times{\bf a}_1)\\
!     {\bf b}_3&=\frac{2\pi}{s}({\bf a}_1\times{\bf a}_2)
!   \end{align*}
!   and finds the unit cell volume $\Omega=|s|$, where
!   $s={\bf a}_1\cdot({\bf a}_2\times{\bf a}_3)$.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
real(8) t1
call r3cross(avec(:,2),avec(:,3),bvec(:,1))
call r3cross(avec(:,3),avec(:,1),bvec(:,2))
call r3cross(avec(:,1),avec(:,2),bvec(:,3))
t1=avec(1,1)*bvec(1,1)+avec(2,1)*bvec(2,1)+avec(3,1)*bvec(3,1)
! compute unit cell volume
omega=abs(t1)
if (omega.lt.1.d-6) then
  write(*,*)
  write(*,'("Error(reciplat) omega too small : ",G18.10)') omega
  write(*,'(" Lattice vectors may be collinear")')
  write(*,*)
  stop
end if
bvec(:,:)=(twopi/t1)*bvec(:,:)
return
end subroutine
!EOC
