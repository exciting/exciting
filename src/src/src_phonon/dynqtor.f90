!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dynqtor (dynq, dynr)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Complex (8), Intent (In) :: dynq (3*natmtot, 3*natmtot, nqpt)
      Complex (8), Intent (Out) :: dynr (3*natmtot, 3*natmtot, &
     & ngridq(1)*ngridq(2)*ngridq(3))
! local variables
      Integer :: ir, iq, i, j, n
      Integer :: isym, lspl, iv (3)
      Integer :: i1, i2, i3, j1, j2, j3
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Real (8) :: s (3, 3), t1
      Complex (8) zt1
! allocatable arrays
      Complex (8), Allocatable :: dyns (:, :)
! external functions
      Real (8) :: r3taxi
      External r3taxi
      Allocate (dyns(3*natmtot, 3*natmtot))
      dynr (:, :, :) = 0.d0
! loop over q-vectors
      Do j1 = 0, ngridq (1) - 1
         v1 (1) = dble (j1) / dble (ngridq(1))
         Do j2 = 0, ngridq (2) - 1
            v1 (2) = dble (j2) / dble (ngridq(2))
            Do j3 = 0, ngridq (3) - 1
               v1 (3) = dble (j3) / dble (ngridq(3))
               iq = iqmap (j1, j2, j3)
! map v1 to the first Brillouin zone
               v2 (:) = v1 (:)
               Call vecfbz (input%structure%epslat, bvec, v2, iv)
! rotate and add the dynamical matrix of the reduced q-point with all symmetries
               n = 0
               dyns (:, :) = 0.d0
               Do isym = 1, nsymcrys
                  lspl = lsplsymc (isym)
                  s (:, :) = dble (symlat(:, :, lspl))
                  Call r3mtv (s, vql(:, iq), v3)
                  Call vecfbz (input%structure%epslat, bvec, v3, iv)
                  If (r3taxi(v2, v3) .Lt. input%structure%epslat) Then
                     Call dynsymapp (isym, vql(:, iq), dynq(:, :, iq), &
                    & dyns)
                     n = n + 1
                  End If
               End Do
               If (n .Eq. 0) Then
                  Write (*,*)
                  Write (*, '("Error(dynqtor): vector ", 3G18.10)') v1
                  Write (*, '(" cannot be mapped to reduced q-point set&
                 &")')
                  Write (*,*)
                  Stop
               End If
               t1 = 1.d0 / dble (n)
               dyns (:, :) = t1 * dyns (:, :)
! loop over R-vectors
               ir = 0
               Do i3 = ngridq (3) / 2 - ngridq (3) + 1, ngridq (3) / 2
                  Do i2 = ngridq (2) / 2 - ngridq (2) + 1, ngridq (2) / &
                 & 2
                     Do i1 = ngridq (1) / 2 - ngridq (1) + 1, ngridq &
                    & (1) / 2
                        ir = ir + 1
                        t1 = twopi * &
                       & (v1(1)*dble(i1)+v1(2)*dble(i2)+v1(3)*dble(i3))
                        zt1 = cmplx (Cos(t1), Sin(t1), 8)
                        Do i = 1, 3 * natmtot
                           Do j = 1, 3 * natmtot
                              dynr (i, j, ir) = dynr (i, j, ir) + zt1 * &
                             & dyns (i, j)
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      t1 = 1.d0 / dble (ngridq(1)*ngridq(2)*ngridq(3))
      dynr (:, :, :) = t1 * dynr (:, :, :)
      Deallocate (dyns)
      Return
End Subroutine
