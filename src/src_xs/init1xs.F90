
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init1xs(vploff)
  use modmain, only: vkloff
  use modxs, only: skipallocs1
  implicit none
  ! arguments
  real(8), intent(in) :: vploff(3)
  ! local variables
  real(8) :: vklofft(3)
  skipallocs1=.true.
  vklofft(:)=vkloff(:)
  vkloff(:)=vploff(:)
  ! call init1 without selected re-allocations and specified k-point offset
  call init1
  skipallocs1=.false.
  vkloff(:)=vklofft(:)
end subroutine init1xs
