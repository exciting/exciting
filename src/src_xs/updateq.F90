
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine updateq(iq)
  use modmain
  use modtddft
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! update the current q-point in the module
  iqcu=iq
  vqlcu(:)=vql(:,iqcu)
end subroutine updateq
