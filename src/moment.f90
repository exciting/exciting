!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: moment
! !INTERFACE:
!
!
Subroutine moment
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total moments by integrating the
!   magnetisation.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, idm
      Real (8) :: sum
! automatic arrays
      Real (8) :: fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax)
      If ( .Not. associated(input%groundstate%spin)) Then
         mommt (:, :) = 0.d0
         mommttot (:) = 0.d0
         momir (:) = 0.d0
         momtot (:) = 0.d0
         Return
      End If
! find the muffin-tin moments
      mommttot (:) = 0.d0
      Do idm = 1, ndmag
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is)
                  fr (ir) = fourpi * magmt (1, ir, ias, idm) * y00 * &
                 & spr (ir, is) ** 2
               End Do
               Call fderiv (-1, nrmt(is), spr(:, is), fr, gr, cf)
               mommt (idm, ias) = gr (nrmt(is))
               mommttot (idm) = mommttot (idm) + mommt (idm, ias)
            End Do
         End Do
      End Do
! find the interstitial moments
      Do idm = 1, ndmag
         sum = 0.d0
         Do ir = 1, ngrtot
            sum = sum + magir (ir, idm) * cfunir (ir)
         End Do
         momir (idm) = sum * omega / dble (ngrtot)
      End Do
      momtot (:) = mommttot (:) + momir (:)
      Return
End Subroutine
!EOC
