!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findkmapkq (vq, voff, map)
      Use modmain
      Use modinput
      Use modxs
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vq (3), voff (3)
      Integer, Intent (Out) :: map (nkpt)
  ! local variables
      Integer :: ik, iv (3), ivt (3)
      Real (8) :: vkq (3)
      Real (8), External :: r3taxi
      Do ik = 1, nkpt
         vkq (:) = vkl (:, ik) + vq (:)
         Call r3frac (input%structure%epslat, vkq, ivt)
         iv (:) = Nint (vkq(:)*input%groundstate%ngridk(:)-voff(:))
         map (ik) = ikmap (iv(1), iv(2), iv(3))
      End Do
End Subroutine findkmapkq
