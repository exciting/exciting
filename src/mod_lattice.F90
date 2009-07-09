

#include "maxdefinitions.inc"
module mod_lattice
!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
real(8)::avec(3, 3)
! inverse of lattice vector matrix
real(8)::ainv(3, 3)
! reciprocal lattice vectors
real(8)::bvec(3, 3)
! inverse of reciprocal lattice vector matrix
real(8)::binv(3, 3)
! unit cell volume
real(8)::omega
! any vector with length less than epslat is considered zero
real(8)::epslat
end module
