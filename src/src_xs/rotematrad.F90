

! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine rotematrad(ngp, igpmap)
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: ngp, igpmap(ngp)
  ! rotate radial integrals
  riaa(:, :, :, :, :, :, :)=riaa(:, :, :, :, :, :, igpmap)
  riloa(:, :, :, :, :, :)=riloa(:, :, :, :, :, igpmap)
  rilolo(:, :, :, :, :)=rilolo(:, :, :, :, igpmap)
end subroutine rotematrad
