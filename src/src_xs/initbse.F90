
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initbse
  use modmain
  use modxs
  implicit none

  ! irreversibly map varialbes specific for BSE (-kernel) to main variables
  nosym=nosymbse
  reducek=reducekbse
  vkloff(:)=vkloffbse(:)
  rgkmax=rgkmaxbse

  ! only one SCF iteration
  maxscl=1

  ! work with regular q-point grid
  qtype='grid'

end subroutine initbse
