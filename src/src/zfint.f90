!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Complex (8) Function zfint (zfmt, zfir)
      Use modmain
      Implicit None
! arguments
      Complex (8), Intent (In) :: zfmt (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zfir (ngrtot)
! local variables
      Integer :: is, ia, ias, ir, irc
      Real (8) :: t1, t2
      Complex (8) zsum
! automatic arrays
      Real (8) :: fr1 (nrcmtmax), fr2 (nrcmtmax), gr (nrcmtmax), cf (3, &
     & nrcmtmax)
      zsum = 0.d0
! interstitial contribution
      Do ir = 1, ngrtot
         zsum = zsum + cfunir (ir) * zfir (ir)
      End Do
      zsum = zsum * omega / dble (ngrtot)
! muffin-tin contribution
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do irc = 1, nrcmt (is)
               t1 = rcmt (irc, is) ** 2
               fr1 (irc) = dble (zfmt(1, irc, ias)) * t1
               fr2 (irc) = aimag (zfmt(1, irc, ias)) * t1
            End Do
            Call fderiv (-1, nrcmt(is), rcmt(:, is), fr1, gr, cf)
            t1 = gr (nrcmt(is))
            Call fderiv (-1, nrcmt(is), rcmt(:, is), fr2, gr, cf)
            t2 = gr (nrcmt(is))
            zsum = zsum + fourpi * y00 * cmplx (t1, t2, 8)
         End Do
      End Do
      zfint = zsum
      Return
End Function
