
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getbsediag
  use modxs
  use m_getunit
  implicit none
  ! local variables
  integer :: un
  real(8) :: re,im
  call getunit(un)
  open(un,file='BSEDIAG.OUT',action='read',form='formatted',status='old')
  read(un,*) re,im
  bsed=cmplx(re,im,8)
  read(un,*) re,im
  bsedl=cmplx(re,im,8)
  read(un,*) re,im
  bsedu=cmplx(re,im,8)
  read(un,*) re,im
  bsedd=cmplx(re,im,8)
  close(un)
end subroutine getbsediag

subroutine getbsedg(fname,iknr,n1,n2,d)
  use modxs
  use m_getunit
  implicit none
  ! arguments
  character(*), intent(in) :: fname
  integer, intent(in) :: iknr,n1,n2
  complex(8), intent(out) :: d(n1,n2)
  ! local variables
  integer :: un,recl,iknrt,n1t,n2t
  call getunit(un)
  inquire(iolength=recl) iknrt,n1t,n2t,d
  open(un,file=trim(fname),action='read',form='unformatted', &
  	status='old',access='direct',recl=recl)
  read(un,rec=iknr) iknrt,n1t,n2t,d
  close(un)
end subroutine getbsedg
