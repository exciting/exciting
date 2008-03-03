
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

character(256) function r2str(r,frmt)
  implicit none
  real(8), intent(in) :: r
  character(*),intent(in)::frmt
  write(r2str,frmt)r
  r2str=trim(adjustl(r2str))
end function r2str
