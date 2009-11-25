!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
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
