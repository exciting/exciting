
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initscr
  use modmain
  use modtddft
  implicit none

  ! irreversibly map varialbes specific for screening to main variables
  nosym=nosymscr
  reducek=reducekscr
!  vkloff(:)=vkloffscr(:)
  rgkmax=rgkmaxscr

  ! only one SCF iteration
  maxscl=1

  ! work with regular q-point grid
  qtype='grid'
  ngridq(:)=ngridk(:)

end subroutine initscr
