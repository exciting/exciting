
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putbsemat(fname,zmat,ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4)
  use m_getunit
  implicit none
  ! arguments
  character(*), intent(in) :: fname
  complex(8), intent(in) :: zmat(n1,n2,n3,n4)
  integer, intent(in) :: ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4
  ! local variables
  integer :: recl,un
  call getunit(un)
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4,zmat
  open(unit=un,file=trim(fname),form='unformatted',action='write', &
       access='direct',recl=recl)
  write(un,rec=ikkp) ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4,zmat
  close(un)
end subroutine putbsemat
