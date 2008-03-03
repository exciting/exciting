
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gettetcw2
  implicit none
contains
  
  subroutine gettetcw2(iq,ik,ib1,ib2,nw,fnam,cw,cwa,cwsurf)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,ib1,ib2,nw
    character(*), intent(in) :: fnam
    real(8), intent(out) :: cw(nw),cwa(nw),cwsurf(nw)
    ! local variables
    integer :: un,irec,recl,iqt
    
    ! q-point
    iqt=iq

    ! record position
    irec=(ik-1)*nstsv*nstsv + (ib1-1)*nstsv + ib2

    ! read from file
    call getunit(un)
    inquire(iolength=recl) cw,cwa,cwsurf
    open(un,file=trim(fnam),form='unformatted',action='read',&
         status='old',access='direct',recl=recl)
    read(un,rec=irec) cw,cwa,cwsurf
    close(un)

  end subroutine gettetcw2

end module m_gettetcw2
