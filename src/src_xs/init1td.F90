
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init1td
  use modxs, only: skipallocs1
  implicit none
  skipallocs1=.true.
  ! call init1 without (re-)allocation of radial functions
  call init1
  skipallocs1=.false.
end subroutine init1td
