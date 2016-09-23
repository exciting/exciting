!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
!
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Character (256) Function stringtim (sec, hrs, d, h, m, s)
      Implicit None
  ! arguments
      Real (8), Intent (In) :: sec, hrs
      Integer, Intent (In) :: d, h, m, s
  ! functions
      Character (256), External :: r2str, i2str
      stringtim = trim (r2str(sec, '(F12.2)')) // ' sec; ' // trim &
     & (r2str(hrs, '(F12.2)')) // ' hrs; ( ' // trim (i2str(d, '(I4)')) &
     & // ' d, ' // trim (i2str(h, '(I3.2)')) // ' h, ' // trim &
     & (i2str(m, '(I3.2)')) // ' m, ' // trim (i2str(s, '(I3.2)')) // '&
     & s )'
End Function stringtim
