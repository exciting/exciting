!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hmlalon (hamilton, is, ia, ngp, apwalm)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (hermiteanMatrix), Intent (Inout) :: hamilton
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
            Do l3 = 0, input%groundstate%lmaxmat
               Do m3 = - l3, l3
                  lm3 = idxlm (l3, m3)
                  Do io = 1, apword (l3, is)
                     zsum = 0.d0
                     Do l2 = 0, input%groundstate%lmaxvr
                        If (Mod(l1+l2+l3, 2) .Eq. 0) Then
                           Do m2 = - l2, l2
                              lm2 = idxlm (l2, m2)
                              zsum = zsum + gntyry (lm1, lm2, lm3) * &
                             & hloa (ilo, io, l3, lm2, ias)
                           End Do
                        End If
                     End Do
                     If (Abs(dble(zsum))+Abs(aimag(zsum)) .Gt. 1.d-20) &
                    & Then
!
! calculate the matrix elements
              !k=((i-1)*i)/2
                        Do j = 1, ngp
                !k=k+1
                           zt1 = zsum * apwalm (j, io, lm3, ias)
!
                           Call Hermiteanmatrix_indexedupdate &
                          & (hamilton, i, j, conjg(zt1))
                        End Do
!
                     End If
                  End Do
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
