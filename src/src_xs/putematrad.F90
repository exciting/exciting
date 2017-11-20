! Copyright(C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine putematrad(iqr, iq)
  use modxs, only: riaa, riloa, rilolo, ematraddir
  use m_genfilname
  use m_getunit
  use modinput

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, iq

  ! Local variables
  character(256) :: fname
  integer :: un

  ! Calculate radial integrals
  call genfilname(basename=trim(adjustl(ematraddir))//'/'//'EMATRAD',&
    & iq=iqr, appfilext=.true., filnam=fname)
  call ematrad(iq)
  call getunit(un)

  open(un, file=trim(fname), form='unformatted', action='write', status='replace')
  write(un) riaa, riloa, rilolo
  close(un)
end subroutine putematrad
