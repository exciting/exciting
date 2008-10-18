
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init1offs(vploff)
  use modmain, only: vkloff
  use modxs, only: init1norealloc
  implicit none
  ! arguments
  real(8), intent(in) :: vploff(3)
  ! local variables
  real(8) :: vt(3)
  logical :: lt
  lt=init1norealloc
  init1norealloc=.true.
  vt(:)=vkloff(:)
  vkloff(:)=vploff(:)
  ! call init1 without selected re-allocations and specified k-point offset
  call init1
  init1norealloc=lt
  vkloff(:)=vt(:)
end subroutine init1offs
