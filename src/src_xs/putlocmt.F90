

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine putlocmt(fname, ik, vk, vq, locmt)
  use modmain
  use m_getunit
  implicit none
  ! arguments
  character(*), intent(in) :: fname
  integer, intent(in) :: ik
  real(8), intent(in) :: vk(3), vq(3)
  complex(8), intent(in) :: locmt(nstfv, nlomax, -lolmax:lolmax, natmtot)
  ! local variables
  integer :: recl, un
  call getunit(un)
  inquire(iolength=recl) vq, vk, nstfv, nlomax, lolmax, locmt
  open(un, file = trim(fname), action = 'write', form = 'unformatted', access = 'direct', &
       recl = recl)
  write(un, rec=ik) vq, vk, nstfv, nlomax, lolmax, locmt
  close(un)
end subroutine putlocmt
