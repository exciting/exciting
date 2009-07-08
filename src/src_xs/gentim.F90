

! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine gentim(sec_in, hrs, days, hours, minutes, seconds)
  implicit none
  ! arguments
  real(8), intent(in) :: sec_in
  real(8), intent(out) :: hrs
  integer, intent(out) :: days, hours, minutes, seconds
  ! local variables
  integer :: isecs
  hrs=sec_in/3600.d0
  isecs=nint(sec_in)
  days=isecs/(60*60*24)
  hours=mod(isecs/(60*60), 24)
  minutes=mod(isecs/60, 60)
  seconds=mod(isecs, 60)
end subroutine gentim
