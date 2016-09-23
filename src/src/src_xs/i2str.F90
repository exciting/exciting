!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Function i2str (i, frmt)
      Implicit None
      Character (256) :: i2str
      Integer, Intent (In) :: i
      Character (*), Intent (In) :: frmt
      Write (i2str, frmt) i
      i2str = trim (adjustl(i2str))
End Function i2str
