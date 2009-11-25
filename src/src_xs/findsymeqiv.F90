!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findsymeqiv (tfbz, vpl, vplr, nsc, sc, ivgsc)
      Use modmain
      Use modinput
      Implicit None
  ! arguments
      Logical, Intent (In) :: tfbz
      Real (8), Intent (In) :: vpl (3), vplr (3)
      Integer, Intent (Out) :: nsc, sc (maxsymcrys), ivgsc (3, &
     & maxsymcrys)
  ! local variables
      Integer :: isym, lspl, iv (3)
      Real (8) :: s (3, 3), v1 (3), t1
      Real (8), External :: r3taxi
  ! symmetries that transform non-reduced q-point to reduced one, namely
  ! vpl = s^-1 * vplr + G_s
      nsc = 0
      Do isym = 1, nsymcrys
         lspl = lsplsymc (isym)
         s (:, :) = dble (symlat(:, :, lspl))
         Call r3mtv (s, vplr, v1)
         If (tfbz) Then
            Call vecfbz (input%structure%epslat, bvec, v1, iv)
         Else
            Call r3frac (input%structure%epslat, v1, iv)
         End If
         t1 = r3taxi (vpl, v1)
         If (t1 .Lt. input%structure%epslat) Then
            nsc = nsc + 1
            sc (nsc) = isym
            ivgsc (:, nsc) = - iv (:)
         End If
      End Do
      If (nsc .Eq. 0) Then
         Write (*,*)
         Write (*, '("Error(findsymeqiv): p-points are not equivalent b&
        &y symmetry")')
         Write (*, '(" vpl  :", 3g18.10)') vpl
         Write (*, '(" vplr :", 3g18.10)') vplr
         Write (*,*)
         Call terminate
      End If
End Subroutine findsymeqiv
