
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine restore1
  use modmain
  use modxs
  implicit none
  nempty=nempty_b
  rgkmax=rgkmax_b
  reducek=reducek
  ngridk(:)=ngridk_b(:)
  vkloff(:)=vkloff_b(:)
  emattype=emattype_b
end subroutine restore1
