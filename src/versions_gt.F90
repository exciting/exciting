
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public
! License. See the file COPYING for license details.

logical function versions_gt(v,vp)
  implicit none
  ! arguments
  integer, intent(in) :: v(3), vp(3)
  versions_gt = &
  (v(1) .gt. vp(1)) .or. &
  ((v(1) .eq. vp(1)) .and. (v(2) .gt. vp(2))) .or. &
  ((v(1) .eq. vp(1)) .and. (v(2) .eq. vp(2)) .and. (v(3) .gt. vp(3)))
end function
