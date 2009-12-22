!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olpalon (overlap, is, ia, ngp, apwalm)
      Use modmain
      Use modfvsystem
      Implicit None
! arguments
      Type (HermiteanMatrix) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
!
!
! local variables
      Integer :: ias, ilo, io, l, m, lm, i, j, k
      Complex (8) zsum
      ias = idxas (ia, is)
      Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         Do m = - l, l
            lm = idxlm (l, m)
            j = ngp + idxlo (lm, ilo, ias)
!
! calculate the matrix elements
      !k=((j-1)*j)/2
            Do i = 1, ngp
       ! k=k+1
               zsum = 0.d0
               Do io = 1, apword (l, is)
                  zsum = zsum + conjg (apwalm(i, io, lm, ias)) * oalo &
                 & (io, ilo, ias)
               End Do
        !o(k)=o(k)+zsum
!
               Call Hermiteanmatrix_indexedupdate (overlap, j, i, zsum)
            End Do
!
         End Do
      End Do
      Return
End Subroutine
