
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xschkstop
  use modxs
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'xschkstop'
  integer :: un
  logical :: exist
  inquire(file='STOP',exist=exist)
  if (exist) then
     call getunit(un)
     write(unitout,'(a)') 'Warning('//thisnam//') stopped with message: '// &
          trim(msg)
     open(un,file='STOP')
     close(un,status='delete')
     call terminate
  end if
end subroutine xschkstop
