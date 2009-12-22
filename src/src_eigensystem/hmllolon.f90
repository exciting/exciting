!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hmllolon (hamilton, is, ia, ngp)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (HermiteanMatrix) :: hamilton
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
               Do m3 = - l3, l3
                  lm3 = idxlm (l3, m3)
                  j = ngp + idxlo (lm3, ilo2, ias)
                  If (i .Le. j) Then
                     zsum = 0.d0
                     Do l2 = 0, input%groundstate%lmaxvr
                        If (Mod(l1+l2+l3, 2) .Eq. 0) Then
                           Do m2 = - l2, l2
                              lm2 = idxlm (l2, m2)
                              zsum = zsum + gntyry (lm1, lm2, lm3) * &
                             & hlolo (ilo1, ilo2, lm2, ias)
                           End Do
                        End If
                     End Do
!
! calculate the matrix elements
!
!
                     Call Hermiteanmatrix_indexedupdate (hamilton, j, &
                    & i, zsum)
!
                  End If
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
