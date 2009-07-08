

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine showunits(un)
  implicit none
  ! arguments
  integer, intent(in) :: un
  ! local variables
  integer :: minu, maxu, u
  logical :: connected
  character(256) :: fname
  minu=-1
  maxu=70000
  write(un, '(a)') 'Info(showunits):'
  write(un, '(a, i9, a, i9, a)') ' opened files between units ', minu, ' and', &
       maxu, ' (unit/name):'
  do u=minu, maxu
     inquire(u, opened=connected)
     if (connected) then
	inquire(u, name=fname)
	write(un, '(i9, 2x, a)') u, trim(fname)
     end if
  end do
end subroutine showunits
