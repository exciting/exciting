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
      Type (HermitianMatrix) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
!
!
! local variables
      Integer :: ias, ilo, io, l, m, lm, i, j, k, lm1, lm2, j1, j2
      Complex (8) zsum
      ias = idxas (ia, is)
      Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         lm1 = idxlm (l,-l)
         lm2 = idxlm (l, l)            
         j1 = ngp + idxlo (lm1, ilo, ias)
         j2 = ngp + idxlo (lm2, ilo, ias)
         Do io = 1, apword (l, is)
           overlap%za(1:ngp,j1:j2)=overlap%za(1:ngp,j1:j2)+conjg(apwalm(:, io, lm1:lm2, ias) * oalo (io, ilo, ias))
         End Do
         do j=j1,j2
            overlap%za(j,1:ngp)=conjg(overlap%za(1:ngp,j))
         End Do
      End Do
      Return
End Subroutine
