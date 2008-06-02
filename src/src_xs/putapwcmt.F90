
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putapwcmt(fname,ik,vk,vq,apwcmt)
  use modmain
  use m_getunit
  implicit none
  ! arguments
  character(*), intent(in) :: fname
  integer, intent(in) :: ik
  real(8), intent(in) :: vk(3),vq(3)
  complex(8), intent(in) :: apwcmt(nstsv,apwordmax,lmmaxapw,natmtot)
  ! local variables
  integer :: recl,un
  call getunit(un)
  inquire(iolength=recl) vq,vk,apwcmt
  open(un,file=trim(fname),action='write',form='unformatted',access='direct', &
       recl=recl)
  write(un,rec=ik) vq,vk,apwcmt
  close(un)
end subroutine putapwcmt
