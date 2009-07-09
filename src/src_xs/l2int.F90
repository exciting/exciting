


! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

integer function l2int(l)
  implicit none
  logical, intent(in) :: l
  if (l) then
     l2int=1
  else
     l2int=0
  end if
end function l2int
