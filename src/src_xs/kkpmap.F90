

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.


subroutine kkpmap(ikkp, nkp, ik, ikp)
  implicit none
  ! arguments
  integer, intent(in) :: ikkp, nkp
  integer, intent(out) :: ik, ikp
  ik=ceiling(0.5d0+nkp-sqrt((0.5d0+nkp)**2-2.d0*ikkp))
  ikp=ikkp+(ik*(ik-1))/2-nkp*(ik-1)
end subroutine kkpmap
