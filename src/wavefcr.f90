!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine wavefcr (lrstp, is, ia, ist, m, ld, wfcr)
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ist
      Integer, Intent (In) :: m
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: wfcr (lmmaxvr, ld, 2)
! local variables
      Integer :: ias, k, l, lm, ir, irc
      Real (8) :: cg1, cg2, t1, t2
      l = spl (ist, is)
      k = spk (ist, is)
      If ((m .Lt.-k) .Or. (m .Gt. k-1)) Go To 10
      ias = idxas (ia, is)
! zero the wavefunction
      irc = 0
      Do ir = 1, nrmt (is), lrstp
         irc = irc + 1
         wfcr (:, irc, :) = 0.d0
      End Do
! calculate the Clebsch-Gordan coefficients
      If (k .Eq. l+1) Then
         cg1 = Sqrt (dble(l+m+1)/dble(2*l+1))
         cg2 = - Sqrt (dble(l-m)/dble(2*l+1))
      Else If (k .Eq. l) Then
         cg1 = Sqrt (dble(l-m)/dble(2*l+1))
         cg2 = Sqrt (dble(l+m+1)/dble(2*l+1))
      Else
         Go To 10
      End If
! detemine the two-component spinor
      irc = 0
      Do ir = 1, nrmt (is), lrstp
         irc = irc + 1
! major component of radial wavefunction
         t1 = rwfcr (ir, 1, ist, ias) / spr (ir, is)
         If (Abs(m) .Le. l) Then
            lm = idxlm (l, m)
            t2 = t1 * cg1
            wfcr (1:lmmaxvr, irc, 1) = t2 * zbshtvr (1:lmmaxvr, lm)
         End If
         If (Abs(m+1) .Le. l) Then
            lm = idxlm (l, m+1)
            t2 = t1 * cg2
            wfcr (1:lmmaxvr, irc, 2) = t2 * zbshtvr (1:lmmaxvr, lm)
         End If
      End Do
      Return
10    Continue
      Write (*,*)
      Write (*, '("Error(wavefcr): mismatched l, k or m : ", 3I4)') l, &
     & k, m
      Write (*, '(" for species ", I4)') is
      Write (*, '(" atom ", I4)') ia
      Write (*, '(" and state ", I6)') ist
      Write (*,*)
      Stop
End Subroutine
