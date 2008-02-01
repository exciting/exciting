
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gettetcw
  implicit none
contains
  
  subroutine gettetcw(iq,ik,iv,ic,nw,fnam,cw,cwa,cwsurf)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,iv,ic,nw
    character(*), intent(in) :: fnam
    real(8), intent(out) :: cw(nw),cwa(nw),cwsurf(nw)
    ! local variables
    integer :: un,irec,recl,iqt
    
    ! q-point
    iqt=iq

    ! record position
    irec=(ik-1)*nstval*nstcon + (iv-1)*nstcon + ic

    ! read from file
    call getunit(un)
    inquire(iolength=recl) cw,cwa,cwsurf
    open(un,file=trim(fnam),form='unformatted',action='read',&
         status='old',access='direct',recl=recl)
    read(un,rec=irec) cw,cwa,cwsurf
    close(un)

  end subroutine gettetcw

end module m_gettetcw
