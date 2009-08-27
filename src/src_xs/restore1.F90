



! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine restore1
  use modmain
use modinput
  use modxs
  implicit none
  input%groundstate%nempty=nempty_b
  input%groundstate%rgkmax=rgkmax_b
  input%groundstate%reducek=reducek_b
  input%groundstate%ngridk(:)=ngridk_b(:)
  input%groundstate%vkloff(:)=vkloff_b(:)
  input%xs%emattype=emattype_b
end subroutine restore1
