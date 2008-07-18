
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getunit
  implicit none
contains
  
  subroutine getunit(un)
    implicit none
    ! parameters
    integer,intent(out) :: un
    ! local variables
    character(*), parameter :: thisnam='getunit'
    integer :: u, u_lo, u_hi
    logical :: connected
    ! lower value for units
    u_lo=100
    ! upper value for units
    u_hi=5000
    do u=u_lo,u_hi
       inquire(u,opened=connected)
       if (.not.connected) then
          un=u
          return
       end if
    end do
    write(*,'("Error(",a,"): no free file unit available between",i6,"and",&
         &i6)') thisnam,u_lo,u_hi
    stop
  end subroutine getunit

end module m_getunit
