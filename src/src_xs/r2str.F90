!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
!
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Character (256) Function r2str (r, frmt)
      Implicit None
      Real (8), Intent (In) :: r
      Character (*), Intent (In) :: frmt
      Write (r2str, frmt) r
      r2str = trim (adjustl(r2str))
End Function r2str
