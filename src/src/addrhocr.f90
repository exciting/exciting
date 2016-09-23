!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: addrhocr
! !INTERFACE:
!
!
Subroutine addrhocr
! !USES:
      Use modmain
! !DESCRIPTION:
!   Adds the core density to the muffin-tin and interstitial densities. A
!   uniform background density is added in the interstitial region to take into
!   account leakage of core charge from the muffin-tin spheres.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir
      Real (8) :: t1, sum1, sum2
! automatic arrays
      Real (8) :: fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax)
      sum1 = 0.d0
      sum2 = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
! add the core density to the muffin-tin density
               rhomt (1, ir, ias) = rhomt (1, ir, ias) + rhocr (ir, &
              & ias) / y00
               fr (ir) = fourpi * rhocr (ir, ias) * spr (ir, is) ** 2
            End Do
! compute the core charge inside the muffin-tins
            Call fderiv (-1, nrmt(is), spr(:, is), fr, gr, cf)
            sum1 = sum1 + gr (nrmt(is))
         End Do
         sum2 = sum2 + dble (natoms(is)) * (4.d0/3.d0) * pi * &
        & (rmt(is)**3)
      End Do
! add remaining core charge to interstitial density
      chgcrlk = chgcr - sum1
      t1 = chgcrlk / (omega-sum2)
      rhoir (:) = rhoir (:) + t1
      Return
End Subroutine
!EOC
