! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine scrwritepmat
  use modxs, only: unitout
  use m_genfilname

  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Calculate momentum matrix elements
  call writepmatxs

  write(unitout, '("Info(scrwritepmat):&
   & momentum matrix elements for screening finished")')
end subroutine scrwritepmat
