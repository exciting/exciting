!
!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rfinp
! !INTERFACE:
!
!
Function rfinp (lrstp, rfmt1, rfmt2, rfir1, rfir2)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rfmt1 : first function in real spherical harmonics for all muffin-tins
!           (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfmt2 : second function in real spherical harmonics for all muffin-tins
!           (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfir1 : first real interstitial function in real-space (in,real(ngrtot))
!   rfir2 : second real interstitial function in real-space (in,real(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two real fuctions over the entire unit cell.
!   The input muffin-tin functions should have angular momentum cut-off
!   {\tt lmaxvr}. In the intersitial region, the integrand is multiplied with
!   the characteristic function, $\tilde{\Theta}({\bf r})$, to remove the
!   contribution from the muffin-tin. See routines {\tt rfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (JKD)
!EOP
!BOC
      Implicit None
      Real (8) :: rfinp
! arguments
      Integer, Intent (In) :: lrstp
      Real (8), Intent (In) :: rfmt1 (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (In) :: rfmt2 (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (In) :: rfir1 (ngrtot)
      Real (8), Intent (In) :: rfir2 (ngrtot)
! local variables
      Integer :: is, ia, ias, ir
      Real (8) :: sum
! external functions
      Real (8) :: rfmtinp
      External rfmtinp
      sum = 0.d0
! interstitial contribution
      Do ir = 1, ngrtot
         sum = sum + rfir1 (ir) * rfir2 (ir) * cfunir (ir)
      End Do
      sum = sum * omega / dble (ngrtot)
! muffin-tin contribution
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            sum = sum + rfmtinp (lrstp, input%groundstate%lmaxvr, &
           & nrmt(is), spr(:, is), lmmaxvr, rfmt1(:, :, ias), rfmt2(:, &
           & :, ias))
         End Do
      End Do
      rfinp = sum
      Return
End Function
!EOC
