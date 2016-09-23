!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findkpt (vpl, isym, ik)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (In) :: vpl (3)
      Integer, Intent (Out) :: isym
      Integer, Intent (Out) :: ik
! local variables
      Integer :: lspl, iv (3)
      Real (8) :: s (3, 3), v1 (3), v2 (3), t1
      Do isym = 1, nsymcrys
         lspl = lsplsymc (isym)
         s (:, :) = dble (symlat(:, :, lspl))
         Call r3mtv (s, vpl, v1)
         Call r3frac (input%structure%epslat, v1, iv)
         Do ik = 1, nkpt
            v2 (:) = vkl (:, ik)
            Call r3frac (input%structure%epslat, v2, iv)
            t1 = Abs (v1(1)-v2(1)) + Abs (v1(2)-v2(2)) + Abs &
           & (v1(3)-v2(3))
            If (t1 .Lt. input%structure%epslat) Return
         End Do
      End Do
      Write (*,*)
      Write (*, '("Error(findkpt): equivalent k-point not in set")')
      Write (*, '(" Requested k-point : ", 3G18.10)') vpl
      Write (*,*)
      Stop
End Subroutine
