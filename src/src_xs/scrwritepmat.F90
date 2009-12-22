!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine scrwritepmat
      Use modxs
      Use m_genfilname
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
  ! calculate momentum matrix elements
      Call writepmatxs
      Write (unitout, '("Info(scrwritepmat): momentum matrix elements f&
     &or screening finished")')
End Subroutine scrwritepmat
