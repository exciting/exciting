!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: symtr2
! !INTERFACE:
!
!
Subroutine symtr2 (t2)
! !USES:
      Use modmain
      Use modxs
! !DESCRIPTION:
!   Symmetrizes a rank-2 tensor with respect to the rotational part of the crystal
!   symmetries:
!   $$ t_{ij}^{\rm sym} = \frac{1}{N_{\alpha}}\sum_{\alpha} \sum_{k,l=1}^3
!     \alpha_{ik} \alpha_{jl} t_{kl}. $$
!   Here, $t_{ij}$ are the components of the rank-2 tensor, $\alpha_{ij}$
!   denotes the rotational part of the crystal symmetries and $N_{\alpha}$ stands
!   for the total number of crystal symmetries.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (Inout) :: t2 (3, 3)
  ! local variables
      Integer :: iop1, iop2, i, j
      Real (8) :: s2 (3, 3)
      s2 (:, :) = 0.d0
      Do iop1 = 1, 3
         Do iop2 = 1, 3
            Do i = 1, 3
               Do j = 1, 3
                  s2 (iop1, iop2) = s2 (iop1, iop2) + symt2 (iop1, &
                 & iop2, i, j) * t2 (i, j)
               End Do
            End Do
         End Do
      End Do
      t2 (:, :) = s2 (:, :)
End Subroutine symtr2
