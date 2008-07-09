
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gettetcw
  implicit none
contains
  
  subroutine gettetcw(iq,ik,i1,i2,n1,n2,nw,fnam,cw,cwa,cwsurf)
    use modmain
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,i1,i2,n1,n2,nw
    character(*), intent(in) :: fnam
    real(8), intent(out) :: cw(nw),cwa(nw),cwsurf(nw)
    ! local variables
    integer :: un,irec,recl,iqt,nstsv_,n1_,n2_
    real(8) :: vql_(3),vkl_(3)
    ! q-point
    iqt=iq
    ! record position
    irec=(ik-1)*n1*n2 + (i1-1)*n2 + i2
    ! read from file
    call getunit(un)
    inquire(iolength=recl) vql_,vkl_,nstsv_,n1_,n2_,cw,cwa,cwsurf
    open(un,file=trim(fnam),form='unformatted',action='read',&
         status='old',access='direct',recl=recl)
    read(un,rec=irec) vql_,vkl_,nstsv_,n1_,n2_,cw,cwa,cwsurf
    close(un)
  end subroutine gettetcw

end module m_gettetcw
