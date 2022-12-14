!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
! Pack/unpack a real function given in the muffin-tin spheres and the 
! interstitial region to a real n-dim array representing n real-space points
Subroutine rfpack (tpack, n, lrstp, rfmt, rfir, nu)
      Use modmain
      Implicit None
! arguments
      Logical, Intent (In) :: tpack
      Integer, Intent (Inout) :: n
      Integer, Intent (In) :: lrstp
      Real (8), Intent (Inout) :: rfmt (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (Inout) :: rfir (ngrtot)
      Real (8), Intent (Out) :: nu (*)
! local variables
      Integer :: is, ia, ias, ir, lm
      If (tpack) Then
! pack the function
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  Do lm = 1, lmmaxvr
                     n = n + 1
                     nu (n) = rfmt (lm, ir, ias)
                  End Do
               End Do
            End Do
         End Do
         Do ir = 1, ngrtot
            n = n + 1
            nu (n) = rfir (ir)
         End Do
      Else
! unpack the function
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  Do lm = 1, lmmaxvr
                     n = n + 1
                     rfmt (lm, ir, ias) = nu (n)
                  End Do
               End Do
            End Do
         End Do
         Do ir = 1, ngrtot
            n = n + 1
            rfir (ir) = nu (n)
         End Do
      End If
      Return
End Subroutine

! Pack/unpack a complex function given in the muffin-tin spheres and the 
! interstitial region to a real 2*n-dim array representing n real-space points
Subroutine zfpack (tpack, n, lmmax, lrstp, nf, zfmt, zfir, nu)
      Use modmain
      Implicit None
! arguments
      Logical, Intent (In) :: tpack
      Integer, Intent (Inout) :: n
      Integer, Intent (In) :: lmmax, lrstp, nf
      Complex (8), Intent (Inout) :: zfmt (lmmax, nrmtmax, natmtot, *)
      Complex (8), Intent (Inout) :: zfir (ngrtot, *)
      Real (8), Intent (Out) :: nu (*)
! local variables
      Integer :: is, ia, ias, ir, lm, jf
      If (tpack) Then
! pack the function
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  Do lm = 1, lmmax
                     Do jf = 1, nf
                        n = n + 2
                        nu (n-1) = dble( zfmt (lm, ir, ias, jf))
                        nu (n) = aimag( zfmt (lm, ir, ias, jf))
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Do ir = 1, ngrtot
            Do jf = 1, nf
               n = n + 2
               nu (n-1) = dble( zfir (ir, jf))
               nu (n) = aimag( zfir (ir, jf))
            End Do
         End Do
      Else
! unpack the function
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  Do lm = 1, lmmax
                     Do jf = 1, nf
                        n = n + 2
                        zfmt (lm, ir, ias, jf) = cmplx( nu (n-1), nu (n), 8)
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Do ir = 1, ngrtot
            Do jf = 1, nf
               n = n + 2
               zfir (ir, jf) = cmplx( nu (n-1), nu (n), 8)
            End Do
         End Do
      End If
      Return
End Subroutine
