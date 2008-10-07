
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

logical function transik(ik,trans)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: ik,trans(3,ndftrans)
  transik=.false.
  ! quick return
  if (trans(1,1).eq.0) then
     transik=.true.
     return
  end if
  if (any(trans(1,:).eq.ik)) transik=.true.
end function transik
