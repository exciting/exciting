!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olplolon (overlap, is, ia, ngp)
      Use modmain
      Use modfvsystem
      Implicit None
! arguments
      Type (hermiteanmatrix), Intent (Inout) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
!
!
! local variables
      Complex (8) :: zt
      Integer :: ias, ilo1, ilo2, l, m, lm, i, j, k
      ias = idxas (ia, is)
      Do ilo1 = 1, nlorb (is)
         l = lorbl (ilo1, is)
         Do ilo2 = 1, nlorb (is)
            If (lorbl(ilo2, is) .Eq. l) Then
               Do m = - l, l
                  lm = idxlm (l, m)
                  i = ngp + idxlo (lm, ilo1, ias)
                  j = ngp + idxlo (lm, ilo2, ias)
                  If (i .Le. j) Then
                     zt = dcmplx (ololo(ilo1, ilo2, ias), 0.0)
                     Call Hermiteanmatrix_indexedupdate (overlap, j, i, &
                    & zt)
                  End If
               End Do
            End If
         End Do
      End Do
      Return
End Subroutine
