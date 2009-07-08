

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine backup1
  use modmain
use modinput
  use modxs
  implicit none
  nempty_b=input%groundstate%nempty
  rgkmax_b=input%groundstate%rgkmax
  reducek_b=input%groundstate%reducek
  ngridk_b(:)=input%groundstate%ngkgrid(:)
  vkloff_b(:)=input%groundstate%vkloff(:)
  emattype_b=input%xs%emattype
end subroutine backup1
