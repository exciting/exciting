
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine updateq(vqlt)
  use modtddft
  implicit none
  ! arguments
  real(8), intent(in) :: vqlt(3)
  ! update the current q-point in the module
  vqlcu(:)=vqlt(:)
end subroutine updateq
