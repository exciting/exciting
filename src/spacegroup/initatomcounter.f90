
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initatomcounters
  use modsymmetries
  use modspacegroup
  implicit none
  integer::is
  nspecies=size( symmetries%WyckoffPositions%wspeciesarray)
  Do is=1, nspecies
    nwpos(is)=size( symmetries%WyckoffPositions%wspeciesarray(is)%wspecies%wposarray)
  end do
end subroutine
