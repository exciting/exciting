!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine zrhogp (gpc, jlgpr, ylmgp, sfacgp, zrhomt, zrhoir, zrho0)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (In) :: gpc
      Real (8), Intent (In) :: jlgpr (0:input%groundstate%lmaxvr, &
     & nrcmtmax, nspecies)
      Complex (8), Intent (In) :: ylmgp (lmmaxvr)
      Complex (8), Intent (In) :: sfacgp (natmtot)
      Complex (8), Intent (In) :: zrhomt (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zrhoir (ngrtot)
      Complex (8), Intent (Out) :: zrho0
! local variables
      Integer :: is, ia, ias
      Integer :: l, m, lm, nrc, ir
      Real (8) :: t1, t2
      Complex (8) zsum1, zsum2
! automatic arrays
      Real (8) :: fr1 (nrcmtmax), fr2 (nrcmtmax)
      Real (8) :: gr (nrcmtmax), cf (3, nrcmtmax)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
! (note that the phase exp(ip.r) is implicit)
      zrho0 = 0.d0
      Do ir = 1, ngrtot
         zrho0 = zrho0 + cfunir (ir) * zrhoir (ir)
      End Do
      zrho0 = zrho0 / dble (ngrtot)
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
! (note that the phase exp(ip.r) is explicit)
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrc
               zsum1 = 0.d0
               lm = 0
               Do l = 0, input%groundstate%lmaxvr
                  lm = lm + 1
                  zsum2 = zrhomt (lm, ir, ias) * ylmgp (lm)
                  Do m = - l + 1, l
                     lm = lm + 1
                     zsum2 = zsum2 + zrhomt (lm, ir, ias) * ylmgp (lm)
                  End Do
                  zsum1 = zsum1 + jlgpr (l, ir, is) * conjg (zil(l)) * &
                 & zsum2
               End Do
               t1 = rcmt (ir, is) ** 2
               fr1 (ir) = dble (zsum1) * t1
               fr2 (ir) = aimag (zsum1) * t1
            End Do
            Call fderiv (-1, nrc, rcmt(:, is), fr1, gr, cf)
            t1 = gr (nrc)
            Call fderiv (-1, nrc, rcmt(:, is), fr2, gr, cf)
            t2 = gr (nrc)
            zrho0 = zrho0 + (fourpi/omega) * conjg (sfacgp(ias)) * &
           & cmplx (t1, t2, 8)
         End Do
      End Do
      Return
End Subroutine
