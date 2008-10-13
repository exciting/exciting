
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getbsediag
  use modxs
  use m_getunit
  implicit none
  ! local variables
  integer :: un
  call getunit(un)
  open(un,file='BSEDIAG.OUT',action='read',form='formatted',status='old')
  read(un,*) bsed
  read(un,*) bsedl
  read(un,*) bsedu
  read(un,*) bsedd
  close(un)
end subroutine getbsediag
