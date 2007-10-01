
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gridsize
! !INTERFACE:
subroutine gridsize
! !USES:
use modmain
! !DESCRIPTION:
!   Finds the ${\bf G}$-vector grid which completely contains the vectors with
!   $G<G_{\rm max}$ and is compatible with the FFT routine. The optimal sizes
!   are given by
!   $$ n_i=\frac{G_{\rm max}|{\bf a}_i|}{\pi}+1, $$
!   where ${\bf a}_i$ is the $i$th lattice vector.
!
! !REVISION HISTORY:
!   Created July 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer i1,i2,i3
real(8) v1(3),v2(3),v3(3),t1,t2
! find optimal grid size for potential and density
ngrid(:)=int(gmaxvr*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
! find next largest FFT-compatible grid size
call nfftifc(ngrid(1))
call nfftifc(ngrid(2))
call nfftifc(ngrid(3))
if ((ngrid(1).le.0).or.(ngrid(2).le.0).or.(ngrid(3).le.0)) then
  write(*,*)
  write(*,'("Error(gridsize): invalid ngrid : ",3I8)') ngrid
  write(*,*)
  stop
end if
! total number of points in grid
ngrtot=ngrid(1)*ngrid(2)*ngrid(3)
! determine integer ranges for grid
intgv(:,1)=ngrid(:)/2-ngrid(:)+1
intgv(:,2)=ngrid(:)/2
! find the total number of G-vectors with |G| < gmaxvr
t1=gmaxvr**2
ngvec=0
do i1=intgv(1,1),intgv(1,2)
  v1(:)=dble(i1)*bvec(:,1)
  do i2=intgv(2,1),intgv(2,2)
    v2(:)=v1(:)+dble(i2)*bvec(:,2)
    do i3=intgv(3,1),intgv(3,2)
      v3(:)=v2(:)+dble(i3)*bvec(:,3)
      t2=v3(1)**2+v3(2)**2+v3(3)**2
      if (t2.lt.t1) ngvec=ngvec+1
    end do
  end do
end do
return
end subroutine
!EOC
