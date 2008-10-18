
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine terminateqry(str)
  use m_getunit
  implicit none
  ! arguments
  character(*) :: str
  ! local variables
  logical :: exis
  integer :: un
  ! check for terminating
  inquire(file='TERMINATE',exist=exis)
  if (exis) then
     write(*,'(a)') 'Error: user termination of program in routine: '// &
          trim(str)
     call getunit(un)
     open(un,file='TERMINATE',action='write')
     close(un,status='delete')
     call terminate
  end if
end subroutine terminateqry
