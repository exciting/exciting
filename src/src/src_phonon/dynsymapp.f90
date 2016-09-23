!
!
!
! Copyright (C) 2005-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dynsymapp (isym, vpl, dyn, dyns)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: isym
      Real (8), Intent (In) :: vpl (3)
      Complex (8), Intent (In) :: dyn (3*natmtot, 3*natmtot)
      Complex (8), Intent (Inout) :: dyns (3*natmtot, 3*natmtot)
! local variables
      Integer :: is, ia, ja, ias, jas
      Integer :: lspl, i, j, k, l, m, n, iv (3)
      Real (8) :: s (3, 3), a (3, 3), b (3, 3), c (3, 3)
      Real (8) :: v1 (3), v2 (3), v3 (3), t1
      Complex (8) zt1
! automatic arrays
      Integer :: map (natmtot)
      Complex (8) zph (natmtot)
! index to spatial rotation in lattice point group
      lspl = lsplsymc (isym)
! check if symmetry is the identity
      If (lspl .Eq. 1) Then
         dyns (:, :) = dyns (:, :) + dyn (:, :)
         Return
      End If
! symmetry in lattice coordinates
      s (:, :) = dble (symlat(:, :, lspl))
! map vpl to the first Brillouin zone
      v1 (:) = vpl (:)
      Call vecfbz (input%structure%epslat, bvec, v1, iv)
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! equivalent atom with this symmetry
            ja = ieqatom (ia, is, isym)
            jas = idxas (ja, is)
            map (ias) = jas
! phase factor
            v2 (:) = input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord(:) + vtlsymc (:, isym)
            Call r3mv (s, v2, v3)
            v3 (:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:) - v3 (:)
            t1 = twopi * (v1(1)*v3(1)+v1(2)*v3(2)+v1(3)*v3(3))
            zph (ias) = cmplx (Cos(t1), Sin(t1), 8)
         End Do
      End Do
! symmetry in Cartesian coordinates
      s (:, :) = symlatc (:, :, lspl)
! rotate and phase-shift dynamical matrix with symmetry
      Do ias = 1, natmtot
         i = 3 * (ias-1)
         k = 3 * (map(ias)-1)
         Do jas = 1, natmtot
            j = 3 * (jas-1)
            l = 3 * (map(jas)-1)
            Do m = 1, 3
               Do n = 1, 3
                  a (m, n) = dble (dyn(i+m, j+n))
                  b (m, n) = aimag (dyn(i+m, j+n))
               End Do
            End Do
            Call r3mtm (s, a, c)
            Call r3mm (c, s, a)
            Call r3mtm (s, b, c)
            Call r3mm (c, s, b)
            zt1 = zph (ias) * conjg (zph(jas))
            Do m = 1, 3
               Do n = 1, 3
                  dyns (k+m, l+n) = dyns (k+m, l+n) + zt1 * cmplx (a(m, &
                 & n), b(m, n), 8)
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
