
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!> Real and reciprocal-space lattice vectors, and volume 
Module mod_lattice
   implicit none 
   !> lattice vectors stored column-wise
   real(8)  :: avec(3,3) ! DN: an alias for input%structure%crystal%basevect (will be initialized in init0
   !> inverse of lattice vector matrix
   Real (8) :: ainv (3, 3)
   !> reciprocal lattice vectors
   Real (8) :: bvec (3, 3)
   !> inverse of reciprocal lattice vector matrix
   Real (8) :: binv (3, 3)
   !> unit cell volume
   Real (8) :: omega
   ! Any vector with length less than epslat is considered zero
   ! replaced by inputstructurereal(8)::epslat
End Module