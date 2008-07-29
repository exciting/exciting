
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initbse
  use modmain
  use modxs
  implicit none
  if (any(vkloffbse.eq.-1.d0)) vkloffbse(:)=vkloff(:)
end subroutine initbse
