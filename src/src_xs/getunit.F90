
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getunit
  implicit none
contains
  
  subroutine getunit(fu)
    implicit none
    ! parameters
    integer,intent(out) :: fu
    ! local variables
    character(*), parameter :: thisnam = 'getunit'
    integer :: u, u_lo, u_hi
    logical :: connected
    u_lo=100
    u_hi=5000
    do u=u_lo,u_hi
       inquire(u,opened=connected)
       if (.not.connected) then
          fu = u
          return
       end if
    end do
    write(*,*) 'Error('//thisnam//'):'
    write(*,*) '  diagnostics: no free file unit found between',u_lo,'and',u_hi
    stop 'Error in getunit'
  end subroutine getunit

end module m_getunit
