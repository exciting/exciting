
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writekmapkq(iq)
  use modmain
  use modxs
  use m_getunit
  implicit none
  integer, intent(in) :: iq
  integer :: un,ik
  call getunit(un)
  open(un,file='KMAPKQ'//trim(filext),form='formatted',action='write', &
       status='replace')
  write(un,'(i9,a)') nkpt, ' : nkpt; k-point, ikmapikq below'
  do ik=1,nkpt
     write(un,'(2i9)') ik,ikmapikq(ik,iq)
  end do
  close(un)
end subroutine writekmapkq
