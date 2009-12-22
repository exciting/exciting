!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genylmgq
! !INTERFACE:
!
!
Subroutine genylmgq (iq, lmax)
! !USES:
      Use modmain
      Use modxs
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
      Implicit None
! arguments
      Integer, Intent (In) :: iq, lmax
! local variables
      Integer :: igq
      Real (8) :: r, tp (2)
      Do igq = 1, ngq (iq)
         Call sphcrd (vgqc(1, igq, iq), r, tp)
         Call genylm (lmax, tp, ylmgq(1, igq, iq))
      End Do
      Return
End Subroutine genylmgq
!EOC
