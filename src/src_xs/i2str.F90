
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

character(256) function i2str(i,frmt)
  implicit none
  integer, intent(in) :: i
  character(*),intent(in)::frmt
  write(i2str,frmt) i
  i2str=trim(adjustl(i2str))
end function i2str
