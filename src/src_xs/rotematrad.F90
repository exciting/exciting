! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rotematrad
! !INTERFACE:
subroutine rotematrad(ngp, igpmap)
! !USES:
  use modxs, only: riaa, riloa, rilolo
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   integer(4), ngp : number of G+q vectors
!   integer(4), igpmap(ngp) : {G+q} --> {G+q_r} index map
!
! !DESCRIPTION:
!   The routine applies the $\left\{ \vec{G}+\vec{q} \right\}$ index
!   map created in {\tt findgqmap} to the radial integral arrays computed for
!   the reduced q set.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich 2016)
!EOP
!BOC

  implicit none

  integer, intent (in) :: ngp, igpmap(ngp)

  ! Rotate radial integrals
  riaa(:, :, :, :, :, :, :) = riaa(:, :, :, :, :, :, igpmap)
  riloa(:, :, :, :, :, :) = riloa(:, :, :, :, :, igpmap)
  rilolo(:, :, :, :, :) = rilolo(:, :, :, :, igpmap)
end subroutine rotematrad
!EOC
