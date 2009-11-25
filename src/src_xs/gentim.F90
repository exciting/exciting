!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gentim (sec_in, hrs, days, hours, minutes, seconds)
      Implicit None
  ! arguments
      Real (8), Intent (In) :: sec_in
      Real (8), Intent (Out) :: hrs
      Integer, Intent (Out) :: days, hours, minutes, seconds
  ! local variables
      Integer :: isecs
      hrs = sec_in / 3600.d0
      isecs = Nint (sec_in)
      days = isecs / (60*60*24)
      hours = Mod (isecs/(60*60), 24)
      minutes = Mod (isecs/60, 60)
      seconds = Mod (isecs, 60)
End Subroutine gentim
