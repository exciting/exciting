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
