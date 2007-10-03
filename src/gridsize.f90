
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
return
end subroutine
!EOC
