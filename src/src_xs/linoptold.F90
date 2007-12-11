
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linoptold
  use modmain
  use modxs
  use modmpi
  use m_getunit
  integer :: un

  if (rank == 0) then
     call linopt
  end if

  call getunit(un)
  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info(linoptold): linear optics (main version) &
       &finished"

end subroutine linoptold
