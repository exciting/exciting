
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phononinterpolate
  use modinput
  implicit none
! initialise universal variables
  Call init0
  Call init2
  call phinterp(input%phonons%interpolate%ngridq, &
    input%phonons%interpolate%vqloff,.true.,.true., &
    input%phonons%interpolate%writeeigenvectors,"PHONON_INTERPOLATE.OUT")
end subroutine
