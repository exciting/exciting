!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine screen
      Use modxs
      Use modinput
      Use m_genfilname
      Implicit None
  ! local variables
      Integer :: nwdft
      nwdft = nwdf
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
  ! call dielectric function with only one frequency point
      Call df
  ! alternative for checking only:
      nwdf = nwdft
      Write (unitout, '(a)') "Info(screen): Screening finished"
End Subroutine screen
