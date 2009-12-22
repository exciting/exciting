!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genylmg
! !INTERFACE:
!
!
Subroutine genylmg
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates a set of spherical harmonics, $Y_{lm}(\hat{\bf G})$, with angular
!   momenta up to {\tt lmaxvr} for the set of ${\bf G}$-vectors.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ig
      Real (8) :: r, tp (2)
! allocate global G-vector spherical harmonic array
      If (allocated(ylmg)) deallocate (ylmg)
      Allocate (ylmg(lmmaxvr, ngvec))
      Do ig = 1, ngvec
         Call sphcrd (vgc(:, ig), r, tp)
         Call genylm (input%groundstate%lmaxvr, tp, ylmg(:, ig))
      End Do
      Return
End Subroutine
!EOC
