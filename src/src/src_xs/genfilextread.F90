!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genfilextread (task)
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: task
      Select Case (task)
      Case (121, 330, 331, 340, 350)
         Call genfilname (iqmt=0, setfilext=.True.)
      Case (430, 440, 441, 445, 450, 451)
         Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
      End Select
End Subroutine genfilextread
