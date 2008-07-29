
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linoptold
  use modmain
  use modxs
  use modmpi
  if (rank.eq.0) then
     call linopt
  end if
  call barrier
  write(unitout,'(a)') "Info(linoptold): linear optics (main version) &
       &finished"
end subroutine linoptold
