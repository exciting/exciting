

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine backup2
  use modmain
  use modxs
  implicit none
  ngridq_b(:)=ngridq(:)
  reduceq_b=reduceq
end subroutine backup2
