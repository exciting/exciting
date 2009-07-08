

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine backup0
  use modmain
use modinput
  use modxs
  implicit none
  filext_b=trim(filext)
  nosym_b=input%groundstate%nosym
end subroutine backup0
