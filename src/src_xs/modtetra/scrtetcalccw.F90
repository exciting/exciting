#ifdef TETRA
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine scrtetcalccw
  use modxs, only: nwdf
  use m_genfilname

  implicit none
  ! Local variables
  integer :: nwdft

  nwdft = nwdf

  ! Only one frequency w=0
  nwdf = 1
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Calculate tetrahedron weights with only one frequency point
  call tetcalccw
  nwdf = nwdft
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

end subroutine scrtetcalccw
!
#endif
