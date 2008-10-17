
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xschkstop
  use modxs
  use m_filedel
  implicit none
  ! local variables
  logical :: exist
  inquire(file='STOP',exist=exist)
  if (exist) then
     write(unitout,'("STOP file exists - stopping with message: ",a)') trim(msg)
     call filedel('STOP')
     call terminate
  end if
end subroutine xschkstop
