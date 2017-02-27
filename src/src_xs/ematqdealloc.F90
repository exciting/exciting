! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematqdealloc
  use modxs, only: evecfv, evecfv0, evalsv0, xih,&
                 & apwcmt, apwcmt0, locmt, locmt0,&
                 & evecfvb, evecfv0b, apwcmtb, apwcmt0b

  implicit none

  deallocate(evecfv, evecfv0, evalsv0)
  deallocate(xih, apwcmt, apwcmt0, locmt, locmt0)
  deallocate(evecfvb, evecfv0b)

end subroutine ematqdealloc
