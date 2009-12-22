!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olpaa (tapp, is, ia, ngp, apwalm, v, o)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Logical, Intent (In) :: tapp
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: v (nmatmax)
      Complex (8), Intent (Inout) :: o (*)
! local variables
      Integer :: ias, l, m, lm, io
      ias = idxas (ia, is)
      Do l = 0, input%groundstate%lmaxmat
         Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
               Call zmatinp (tapp, ngp, zhalf, apwalm(:, io, lm, ias), &
              & apwalm(:, io, lm, ias), v, o)
            End Do
         End Do
      End Do
      Return
End Subroutine
