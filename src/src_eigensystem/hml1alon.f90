!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hml1alon (hamilton, is, ia, ngp, apwalm)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (hermitianMatrix), Intent (Inout) :: hamilton
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
!
!
! local variables
      Integer :: ias, io, ilo, i, j, k
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Complex (8) zsum, zt1
      ias = idxas (ia, is)
      Do ilo = 1, nlorb (is)
         l1 = lorbl (ilo, is)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            i = ngp + idxlo (lm1, ilo, ias)
            if (lm1.le.input%groundstate%lmaxmat) then
              Do io = 1, apword (l1, is)
                zsum = gntyry (lm1, 1, lm1) * h1loa (ilo, io, l1, ias)
                If (Abs(dble(zsum))+Abs(aimag(zsum)) .Gt. 1.d-20) Then
                  Do j = 1, ngp
                    zt1 = zsum * apwalm (j, io, lm1, ias)
                    Call Hermitianmatrix_indexedupdate (hamilton, i, j, conjg(zt1))
                  EndDo
                endif
              enddo
            endif
         End Do
      End Do
      Return
End Subroutine
