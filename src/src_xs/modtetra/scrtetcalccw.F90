#ifdef TETRA
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine scrtetcalccw
      Use modxs
      Use m_genfilname
      Implicit None
  ! local variables
      Integer :: nwdft
      nwdft = nwdf
  ! only one frequency w=0
      nwdf = 1
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
  ! calculate tetrahedron weights with only one frequency point
      Call tetcalccw
      nwdf = nwdft
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
End Subroutine scrtetcalccw
!
#endif
