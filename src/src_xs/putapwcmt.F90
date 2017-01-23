! Copyright(C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine putapwcmt(fname, ik, vk, vq, apwcmt)
  use modmain
  use modinput
  use m_getunit

  implicit none

  ! Arguments
  character(*), intent(in) :: fname
  integer, intent(in) :: ik
  real(8), intent(in) :: vk(3), vq(3)
  complex(8), intent(in) :: apwcmt(nstfv, apwordmax, lmmaxapw, natmtot)

  ! local variables
  integer :: reclen, un

  call getunit(un)

  inquire(iolength=reclen) vq, vk, nstfv, apwordmax,&
    & input%groundstate%lmaxapw, apwcmt

  open(un, file=trim(fname), action='write', form='unformatted',&
    & access='direct', recl=reclen)

  write(un, rec=ik) vq, vk, nstfv, apwordmax, input%groundstate%lmaxapw, apwcmt

  close(un)
end subroutine putapwcmt
