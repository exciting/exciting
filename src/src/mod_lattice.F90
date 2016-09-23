
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_lattice
!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
!replaced by inputstructurereal(8)::avec(3, 3)
! inverse of lattice vector matrix
      Real (8) :: ainv (3, 3)
! reciprocal lattice vectors
      Real (8) :: bvec (3, 3)
! inverse of reciprocal lattice vector matrix
      Real (8) :: binv (3, 3)
! unit cell volume
      Real (8) :: omega
! any vector with length less than epslat is considered zero
!replaced by inputstructurereal(8)::epslat
End Module
!
