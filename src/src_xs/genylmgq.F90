


! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genylmgq
! !INTERFACE:


subroutine genylmgq(iq, lmax)
! !USES:
use modmain
use modxs
! !DESCRIPTION:
!   Generates a set of spherical harmonics, $Y_{lm}(\widehat{{\bf G}+{\bf q}})$,
!   with angular
!   momenta up to {\tt lmax} for the set of ${\bf G+q}$-vectors. Based upon
!   the routine genylmg.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: iq, lmax
! local variables
integer::igq
real(8)::r, tp(2)
do igq=1, ngq(iq)
  call sphcrd(vgqc(1, igq, iq), r, tp)
  call genylm(lmax, tp, ylmgq(1, igq, iq))
end do
return
end subroutine genylmgq
!EOC
