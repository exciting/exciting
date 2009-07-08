
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine backup0
  use modmain
  use modxs
  implicit none
  filext_b=trim(filext)
  nosym_b=nosym
  swidth_b=swidth
  lmaxapw_b=lmaxapw
  lmaxmat_b=lmaxmat
end subroutine backup0
