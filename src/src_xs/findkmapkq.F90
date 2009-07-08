

! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine findkmapkq(vq, voff, map)
  use modmain
use modinput
  use modxs
  implicit none
  ! arguments
  real(8), intent(in) :: vq(3), voff(3)
  integer, intent(out) :: map(nkpt)
  ! local variables
  integer :: ik, iv(3), ivt(3)
  real(8) :: vkq(3)
  real(8), external :: r3taxi
  do ik=1, nkpt
     vkq(:)=vkl(:, ik)+vq(:)
     call r3frac(input%structure%epslat, vkq, ivt)
     iv(:)=nint(vkq(:)*input%groundstate%ngkgrid(:)-voff(:))
     map(ik)=ikmap(iv(1), iv(2), iv(3))
  end do
end subroutine findkmapkq
