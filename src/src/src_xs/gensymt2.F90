!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gensymt2 (maxsymcrys, nsymcrys, symlatc, lsplsymc, symt2)
      Implicit None
  ! arguments
      Integer, Intent (In) :: maxsymcrys
      Integer, Intent (In) :: nsymcrys
      Real (8) :: symlatc (3, 3, 48)
      Integer, Intent (In) :: lsplsymc (maxsymcrys)
      Real (8), Intent (Out) :: symt2 (3, 3, 3, 3)
  ! local variables
      Integer :: isym, i, j, iop1, iop2
      Real (8) :: s (3, 3), sc (3, 3)
      Do iop1 = 1, 3
         Do iop2 = 1, 3
            s (:, :) = 0.d0
            Do isym = 1, nsymcrys
               sc (:, :) = dble (symlatc(:, :, lsplsymc(isym)))
               Do i = 1, 3
                  Do j = 1, 3
                     s (i, j) = s (i, j) + sc (i, iop1) * sc (j, iop2)
                  End Do
               End Do
            End Do
            symt2 (iop1, iop2, :, :) = s (:, :) / dble (nsymcrys)
         End Do
      End Do
End Subroutine gensymt2
