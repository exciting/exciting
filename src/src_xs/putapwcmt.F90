! Copyright(C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine putapwcmt(fname, ik, vk, vq, apwcmt)
  use modmain
  use modinput
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use precision, only: i32, long_int, dp

  implicit none

  ! Arguments
  character(*), intent(in) :: fname
  integer(i32), intent(in) :: ik
  real(dp), intent(in) :: vk(3), vq(3)
  complex(dp), intent(in) :: apwcmt(nstfv, apwordmax, lmmaxapw, natmtot)

  ! local variables
  integer(i32) :: un
  integer(long_int) :: reclen

  call inquire_large( reclen, vq, vk, [nstfv, apwordmax, input%groundstate%lmaxapw], apwcmt )
  call open_direct_unformatted_large( un, trim( adjustl( fname ) ), "write", reclen, "unknown" )

  write(un, rec=ik) vq, vk, nstfv, apwordmax, input%groundstate%lmaxapw, apwcmt

  close(un)
end subroutine putapwcmt
