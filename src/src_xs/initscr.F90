
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initscr
  use modmain
  use modxs
  implicit none
  if (any(ngridkscr.eq.-1)) ngridkscr(:)=ngridk(:)
  if (any(vkloffscr.eq.-1.d0)) vkloffscr(:)=vkloff(:)
  if (nemptyscr.eq.-1) nemptyscr=nempty
  if (rgkmaxscr.eq.-1.d0) rgkmaxscr=rgkmax
end subroutine initscr
