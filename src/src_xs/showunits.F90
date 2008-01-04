
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine showunits
  ! *** Debug routine ***
  implicit none
  integer :: minu,maxu,u
  logical :: connected
  character(256) :: fname
  minu=-1
  maxu=5000
  write(*,'(a)') 'Info(showunits):'
  write(*,'(a,i9,a,i9,a)') ' opened files between units ',minu,' and', &
       maxu,' (unit/name):'
  do u=minu,maxu
     inquire(u,opened=connected)
     if (connected) then
        inquire(u,name=fname)
        write(*,*) u,trim(fname)
     end if
  end do
end subroutine showunits
