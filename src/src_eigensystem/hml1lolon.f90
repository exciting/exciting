!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hml1lolon (hamilton, is, ia, ngp)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (HermitianMatrix) :: hamilton
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
!
!
! local variables
      Integer :: ias, ilo1, ilo2, i, j, k
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Complex (8) zsum
      ias = idxas (ia, is)
      Do ilo1 = 1, nlorb (is)
         l1 = lorbl (ilo1, is)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            i = ngp + idxlo (lm1, ilo1, ias)
            Do ilo2 = 1, nlorb (is)
               l3 = lorbl (ilo2, is)
               if (l3.eq.l1) then
                  j = ngp + idxlo (lm1, ilo2, ias)
                  If (i .Le. j) Then
                    zsum = gntyry (lm1, 1, lm1) * h1lolo (ilo1, ilo2, ias)
                    Call Hermitianmatrix_indexedupdate (hamilton, j, i, zsum)
                  End If
	       endif
            End Do
         End Do
      End Do
      Return
End Subroutine
