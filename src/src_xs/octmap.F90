
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

integer function octmap(i1,i2)
  implicit none
  integer, intent(in) :: i1,i2
  ! second index runs fastest
  octmap=3*(i1-1)+i2
end function octmap
