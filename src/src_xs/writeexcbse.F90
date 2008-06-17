
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeexcbse(iq,n,w,eps,fn)
  use modmain
  use modxs
  use m_getunit
  use m_writevars
  implicit none
  ! arguments
  integer, intent(in) :: iq,n
  real(8), intent(in) :: w(n)
  complex(8), intent(in) :: eps(n)
  character(*), intent(in) :: fn
  ! local variables
  integer :: iw,un
  call getunit(un)
  open(un,file=trim(fn),action='write')
!  write(un,'(4g18.10)') (w(iw)*escale,eps(iw),kkeps(iw),iw=1,n)
  ! write relevant parameters to file
  call writevars(un,iq)
  close(un)
end subroutine writeexcbse
