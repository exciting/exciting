!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: findsymi
! !INTERFACE:
!
!
Subroutine findsymi (epslat, maxsymcrys, nsymcrys, symlat, lsplsymc, &
& vtlsymc, isymlat, scimap)
! !USES:
  use modmpi
! !DESCRIPTION:
!   Throughout the code the symmetries are understood to be applied in a way
!   $$ (\alpha_S|\alpha_R|{\bf t}) {\bf x} = \alpha_S\alpha_R
!   ({\bf x}+{\bf t})$$
!   which is different from the commonly used definition
!   $\{\alpha|\tau\}x=\alpha x+\tau$ -- see routine {\tt findsymcrys}.
!   This difference affects the inverse of the fractional translation
!   but has no effect on the inverse of the rotational part, so the inverse
!   spacegroup symmetry operations are the same for both definitions.
!
! !REVISION HISTORY:
!   Created April 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (In) :: epslat
      Integer, Intent (In) :: symlat (3, 3, 48)
      Integer, Intent (In) :: maxsymcrys, nsymcrys
      Integer, Intent (In) :: lsplsymc (nsymcrys)
      Real (8), Intent (In) :: vtlsymc (3, maxsymcrys)
      Integer, Intent (In) :: isymlat (48)
      Integer, Intent (Out) :: scimap (maxsymcrys)
  ! local variables
      Real (8) :: c (3, 3), si (3, 3), sj (3, 3), vtl (3)
      Integer :: i, isym, jsym, lspli, lsplj, iv (3)
      scimap (:) = 0
      Do isym = 1, nsymcrys
         lspli = lsplsymc (isym)
         si (:, :) = dble (symlat(:, :, lspli))
         Do jsym = 1, nsymcrys
            lsplj = lsplsymc (jsym)
            sj (:, :) = dble (symlat(:, :, lsplj))
        ! translation
            vtl (:) = vtlsymc (:, jsym)
            vtl = matmul (sj, vtl)
            vtl (:) = vtl (:) + vtlsymc (:, isym)
            Call r3frac (epslat, vtl, iv)
        ! rotation
            Call r3mm (si, sj, c)
        ! subract unit matrix
            Forall (i=1:3)
               c (i, i) = c (i, i) - 1.d0
            End Forall
            If ((sum(vtl) .Lt. epslat) .And. (sum(Abs(c)) .Lt. epslat)) &
           & Then
           ! isym is inverse of jsym
               scimap (isym) = jsym
               Go To 10
            End If
         End Do
10       Continue
     ! check if inverse symmetry is consistent with inverse lattice symmetry
         If (isymlat(lsplsymc(isym)) .Ne. lsplsymc(jsym)) Then
            Write (*,*)
            Write (*, '("Error(findsymi): inconsistency with inverse la&
           &ttice symmetry")')
            Write (*, '(" space group symmetry:", t40, i6)') isym
            Write (*, '(" lattice symmetry:", t40, i6)') lsplsymc &
           & (isym)
            Write (*, '(" inverse lattice symmetry:", t40, i6)') &
           & isymlat (lsplsymc(isym))
            Write (*, '(" proposed inverse lattice symmetry:", t40, i6)&
           &') lsplsymc (jsym)
            Write (*,*)
            Call terminate
         End If
      End Do
End Subroutine findsymi
!EOC
