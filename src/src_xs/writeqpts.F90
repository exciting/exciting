! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writeqpts
! !INTERFACE:
!
subroutine writeqpts
! !USES:
  use mod_qpoint, only: nqpt, vql, vqc
  use mod_misc, only: task
  use modxs, only: ngq, nqptr, vqlr, vqcr, wqptr
  use m_getunit
  use m_genfilname
! !DESCRIPTION:
!   Writes the ${\bf q}$-points in lattice coordinates, weights and number of
!   ${\bf G+q}$-vectors to the file {\tt QPOINTS.OUT}. Based on the routine
!   {\tt writekpts}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC

  implicit none

  ! Local variables
  integer :: iq, un
  character(256) :: filnam

  call getunit(un)
  Call genfilname(basename='QPOINTS', appfilext=.True., filnam=filnam)

  open(un, file=trim(filnam), action='WRITE', form='FORMATTED')
  write(un, '(I6, " : nqpt; q-point, vql, vqc, ngq below")') nqpt

  do iq = 1, nqpt
    write(un, '(i6, 6E18.10, i8)') iq, vql(:, iq), vqc(:, iq), ngq(iq)
  end do

  close(un)

  ! Write out reduced q-point set for screened Coulomb interaction
  if(task .eq. 440) then
    call genfilname(basename='QPOINTSR', appfilext=.True., filnam=filnam)

    open(un, file=trim(filnam), action='write', form='formatted', status='replace')
    write(un, '(i6, " : nqptr; q-point, vqlr, vqcr, wqptr below")') nqptr

    do iq = 1, nqptr
      write(un, '(i6, 7E18.10)') iq, vqlr(:, iq), vqcr(:, iq), wqptr(iq)
    end do

    close(un)

  end if
end subroutine writeqpts
!EOC
