
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine symt2app (oct1, oct2, n, symt2, t, tsym)
      Implicit None
  ! arguments
      Integer, Intent (In) :: oct1, oct2, n
      Real (8), Intent (In) :: symt2 (3, 3, 3, 3)
      Complex (8), Intent (In) :: t (3, 3, n)
      Complex (8), Intent (Out) :: tsym (n)
  ! local variables
      Integer :: i, j
  ! symmetrize the macroscopic dielectric function tensor
      tsym (:) = (0.d0, 0.d0)
      Do i = 1, 3
         Do j = 1, 3
            tsym (:) = tsym (:) + symt2 (oct1, oct2, i, j) * t (i, j, :)
         End Do
      End Do
End Subroutine
