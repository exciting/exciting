

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine chkptchk
  use modxs
  implicit none
  ! local variables
  logical :: exis
  inquire(file=trim(fnresume), exist=exis)
  if (exis) then
     write(*, *)
     write(*, '("Error(chkptchk): stale checkpoint file found.")')
     write(*, '(" Either your previous calculation crashed or another")')
     write(*, '(" instance of exciting is already running")')
     write(*, *)
     call terminate
  end if
end subroutine chkptchk
