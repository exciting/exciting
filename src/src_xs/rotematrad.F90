!
!
!
! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rotematrad (ngp, igpmap)
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: ngp, igpmap (ngp)
  ! rotate radial integrals
      riaa (:, :, :, :, :, :, :) = riaa (:, :, :, :, :, :, igpmap)
      riloa (:, :, :, :, :, :) = riloa (:, :, :, :, :, igpmap)
      rilolo (:, :, :, :, :) = rilolo (:, :, :, :, igpmap)
End Subroutine rotematrad
