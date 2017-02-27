! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Subroutine findkmapkq(vq, voff, map)
  Use modinput
  use mod_kpoint, only: nkpt, vkl, ikmap

  Implicit None

  Real(8), Intent(In) :: vq(3), voff(3)
  Integer, Intent(Out) :: map(nkpt)

  ! local variables
  Integer :: ik, iv(3), ivt(3)
  Real(8) :: vkq(3)

  Do ik = 1, nkpt

     ! k' = k+q
     vkq(:) = vkl(:, ik) + vq(:)

     ! Map k' back to [0,1)
     Call r3frac(input%structure%epslat, vkq, ivt)

     ! Get 3d index of k'
     iv(:) = Nint(vkq(:)*input%groundstate%ngridk(:)-voff(:))

     ! Map 3d index to 1d index
     map(ik) = ikmap(iv(1), iv(2), iv(3))

  End Do
End Subroutine findkmapkq
