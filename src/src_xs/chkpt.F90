
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine chkpt(ncpt,cptv,mesg)
  use modxs
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: ncpt,cptv(ncpt)
  character(*), intent(in) :: mesg
  ! local variables
  integer :: un
  character(256) :: str
  write(str,*) ncpt
  call getunit(un)
  open(un,file=trim(fnresume),form='formatted',action='write',status='replace')
  write(un,'(i8," : length of checkpoint vector")') ncpt
  write(un,'('//trim(adjustl(str))//'i8," : checkpoint vector")') cptv(:)
  write(un,'(" (",a,")")') trim(mesg)
  close(un)
end subroutine chkpt
