!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olplolo (tapp, is, ia, ngp, v, o)
      Use modmain
      Implicit None
! arguments
      Logical, Intent (In) :: tapp
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: v (nmatmax)
      Complex (8), Intent (Inout) :: o (*)
! local variables
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
                     If (tapp) Then
! apply the overlap operator to v
                        o (i) = o (i) + ololo (ilo1, ilo2, ias) * v (j)
                        If (i .Ne. j) o (j) = o (j) + ololo (ilo1, &
                       & ilo2, ias) * v (i)
                     Else
! calculate the matrix elements
                        k = i + ((j-1)*j) / 2
                        o (k) = o (k) + ololo (ilo1, ilo2, ias)
                     End If
                  End If
               End Do
            End If
         End Do
      End Do
      Return
End Subroutine
