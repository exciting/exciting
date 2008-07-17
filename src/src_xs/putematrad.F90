
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putematrad(iq)
  use modmain
  use modxs
  use m_genfilname
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  character(256) :: fname
  integer :: un
  ! calculate radial integrals
  call genfilname(basename='EMATRAD',iq=iq,filnam=fname)
  call ematrad(iq )
  call getunit(un)
  open(un,file=trim(fname),form='unformatted',action='write', &
       status='replace')
  write(un) riaa,riloa,rilolo
  close(un)
end subroutine putematrad
