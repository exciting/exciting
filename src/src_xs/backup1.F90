
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine backup1
  use modmain
  use modxs
  implicit none
  nempty_b=nempty
  rgkmax_b=rgkmax
  reducek_b=reducek
  ngridk_b(:)=ngridk(:)
  vkloff_b(:)=vkloff(:)
end subroutine backup1
