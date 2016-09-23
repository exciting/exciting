!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Real (8) Function rfint (rfmt, rfir)
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: rfmt (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (In) :: rfir (ngrtot)
! local variables
      Integer :: is, ia, ias, ir
      Real (8) :: sum
! automatic arrays
      Real (8) :: fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax)
      sum = 0.d0
! interstitial contribution
      Do ir = 1, ngrtot
         sum = sum + rfir (ir) * cfunir (ir)
      End Do
      sum = sum * omega / dble (ngrtot)
! muffin-tin contribution
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               fr (ir) = rfmt (1, ir, ias) * spr (ir, is) ** 2
            End Do
            Call fderiv (-1, nrmt(is), spr(:, is), fr, gr, cf)
            sum = sum + gr (nrmt(is))
         End Do
      End Do
      rfint = sum
      Return
End Function
