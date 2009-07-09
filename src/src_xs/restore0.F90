
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine restore0
  use modmain
  use modxs
  implicit none
  filext=trim(filext_b)
  nosym=nosym_b
  swidth=swidth_b
  lmaxapw=lmaxapw_b
  lmaxmat=lmaxmat_b
end subroutine restore0
