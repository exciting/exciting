!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rvfcross
! !INTERFACE:
!
!
Subroutine rvfcross (rvfmt1, rvfmt2, rvfir1, rvfir2, rvfmt3, rvfir3)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   rvfmt1 : first input muffin-tin field (in,real(lmmaxvr,nrmtmax,natmtot,3))
!   rvfmt2 : second input muffin-tin field (in,real(lmmaxvr,nrmtmax,natmtot,3))
!   rvfir1 : first input interstitial field (in,real(ngrtot,3))
!   rvfir2 : second input interstitial field (in,real(ngrtot,3))
!   rvfmt3 : output muffin-tin field (out,real(lmmaxvr,nrmtmax,natmtot,3))
!   rvfir3 : output interstitial field (out,real(ngrtot,3))
! !DESCRIPTION:
!   Given two real vector fields, ${\bf f}_1$ and ${\bf f}_2$, defined over the
!   entire unit cell, this routine computes the local cross product
!   $$ {\bf f}_3({\bf r})\equiv{\bf f}_1({\bf r})\times{\bf f}_2({\bf r}). $$
!
! !REVISION HISTORY:
!   Created February 2007 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: rvfmt1 (lmmaxvr, nrmtmax, natmtot, 3)
      Real (8), Intent (In) :: rvfmt2 (lmmaxvr, nrmtmax, natmtot, 3)
      Real (8), Intent (In) :: rvfir1 (ngrtot, 3)
      Real (8), Intent (In) :: rvfir2 (ngrtot, 3)
      Real (8), Intent (Out) :: rvfmt3 (lmmaxvr, nrmtmax, natmtot, 3)
      Real (8), Intent (Out) :: rvfir3 (ngrtot, 3)
! local variables
      Integer :: is, ia, ias, ir, itp, i
      Real (8) :: v1 (3), v2 (3), v3 (3)
! automatic arrays
      Real (8) :: rftp1 (lmmaxvr, 3), rftp2 (lmmaxvr, 3)
! muffin-tin region
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               Do i = 1, 3
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 & lmmaxvr, rvfmt1(:, ir, ias, i), 1, 0.d0, rftp1(:, &
                 & i), 1)
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 & lmmaxvr, rvfmt2(:, ir, ias, i), 1, 0.d0, rftp2(:, &
                 & i), 1)
               End Do
               Do itp = 1, lmmaxvr
                  v1 (:) = rftp1 (itp, :)
                  v2 (:) = rftp2 (itp, :)
                  Call r3cross (v1, v2, v3)
                  rftp1 (itp, :) = v3 (:)
               End Do
               Do i = 1, 3
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                 & lmmaxvr, rftp1(:, i), 1, 0.d0, rvfmt3(:, ir, ias, &
                 & i), 1)
               End Do
            End Do
         End Do
      End Do
! interstitial region
      Do ir = 1, ngrtot
         v1 (:) = rvfir1 (ir, :)
         v2 (:) = rvfir2 (ir, :)
         Call r3cross (v1, v2, v3)
         rvfir3 (ir, :) = v3 (:)
      End Do
      Return
End Subroutine
!EOC
