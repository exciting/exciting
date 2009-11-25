!
!
!
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqdealloc
      Use modmain
      Use modmpi
      Use modxs
      Implicit None
      Deallocate (evecfv, evecfv0, evalsv0, xih, apwcmt, apwcmt0, &
     & locmt, locmt0)
End Subroutine ematqdealloc
