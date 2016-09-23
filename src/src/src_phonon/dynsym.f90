!
!
!
! Copyright (C) 2006-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dynsym (vpl, dynp)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (In) :: vpl (3)
      Complex (8), Intent (Inout) :: dynp (3*natmtot, 3*natmtot)
! local variables
      Integer :: iv (3), isym, lspl, i, j, n
      Real (8) :: v1 (3), v2 (3), s (3, 3), t1
! automatic arrays
      Complex (8) dyns (3*natmtot, 3*natmtot)
! external functions
      Real (8) :: r3taxi
      External r3taxi
! map input vector to first Brillouin zone
      v1 (:) = vpl (:)
      Call vecfbz (input%structure%epslat, bvec, v1, iv)
      n = 0
      dyns (:, :) = 0.d0
! use the symmetries which leave vpl invariant
      Do isym = 1, nsymcrys
         lspl = lsplsymc (isym)
         s (:, :) = dble (symlat(:, :, lspl))
         Call r3mtv (s, v1, v2)
         Call vecfbz (input%structure%epslat, bvec, v2, iv)
         If (r3taxi(v1, v2) .Lt. input%structure%epslat) Then
            Call dynsymapp (isym, v1, dynp, dyns)
            n = n + 1
         End If
      End Do
      If (n .Eq. 0) Then
         Write (*,*)
         Write (*, '("Error(dynsym): no symmetries leave vpl invariant"&
        &)')
         Write (*,*)
         Stop
      End If
      t1 = 1.d0 / dble (n)
      dynp (:, :) = t1 * dyns (:, :)
! make the matrix Hermitian
      Do i = 1, 3 * natmtot
         Do j = i, 3 * natmtot
            dynp (i, j) = 0.5d0 * (dynp(i, j)+conjg(dynp(j, i)))
            dynp (j, i) = conjg (dynp(i, j))
         End Do
      End Do
      Return
End Subroutine
