
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initbse
  use modmain
  use modtddft
  implicit none

!!$  ! irreversibly map varialbes specific for BSE (-kernel) to main variables
!!$  nosym=nosymbse
!!$  reducek=reducekbse
!!$  ngridk(:)=ngridkbse(:)
!!$  vkloff(:)=vkloffbse(:)
!!$  rgkmax=rgkmaxbse
!!$  nempty=nemptybse
!!$
!!$  ! only one SCF iteration
!!$  maxscl=1
!!$
!!$  ! work with regular q-point grid
!!$  qtype='grid'

  ! initialize number of empty states
  nempty=nstuoccbse-1

end subroutine initbse
