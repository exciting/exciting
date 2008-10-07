
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genylmg
! !INTERFACE:
subroutine genylmg
! !USES:
use modmain
! !DESCRIPTION:
!   Generates a set of spherical harmonics, $Y_{lm}(\hat{\bf G})$, with angular
!   momenta up to {\tt lmaxvr} for the set of ${\bf G}$-vectors.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ig
real(8) r,tp(2)
! allocate global G-vector spherical harmonic array
if (allocated(ylmg)) deallocate(ylmg)
allocate(ylmg(lmmaxvr,ngvec))
do ig=1,ngvec
  call sphcrd(vgc(:,ig),r,tp)
  call genylm(lmaxvr,tp,ylmg(:,ig))
end do
return
end subroutine
!EOC

