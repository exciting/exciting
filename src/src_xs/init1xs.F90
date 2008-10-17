
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init1xs(vploff)
  use modmain, only: vkloff
  use modxs, only: init1norealloc
  implicit none
  ! arguments
  real(8), intent(in) :: vploff(3)
  ! local variables
  real(8) :: vklofft(3)
  init1norealloc=.true.
  vklofft(:)=vkloff(:)
  vkloff(:)=vploff(:)
  ! call init1 without selected re-allocations and specified k-point offset
  call init1
  init1norealloc=.false.
  vkloff(:)=vklofft(:)
end subroutine init1xs
