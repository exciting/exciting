

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine restore0
  use modmain
use modinput
  use modxs
  implicit none
  filext=trim(filext_b)
  input%groundstate%nosym=nosym_b
end subroutine restore0
