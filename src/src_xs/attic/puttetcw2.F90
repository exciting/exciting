
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_puttetcw2
  implicit none
contains

  subroutine puttetcw2(iq,ik,ib1,ib2,filnam,cw,cwa,cwsurf)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,ib1,ib2
    character(*), intent(in) :: filnam
    real(8), intent(in) :: cw(:),cwa(:),cwsurf(:)
    ! local variables
    integer :: un, recl, irec, iqt

    ! q-point
    iqt=iq

    ! record position
    irec=(ik-1)*nstsv*nstsv + (ib1-1)*nstsv + ib2

    ! I/O record length
    inquire(iolength=recl) cw,cwa,cwsurf
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='write',access='direct',recl=recl)
    write(un,rec=irec) cw,cwa,cwsurf
    close(un)

  end subroutine puttetcw2

end module m_puttetcw2
