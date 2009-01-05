#include "maxdefinitions.inc"
module mod_phonon
!--------------------------!
!     phonon variables     !
!--------------------------!
! number of primitive unit cells in phonon supercell
integer nphcell
! Cartesian offset vectors for each primitive cell in the supercell
real(8) vphcell(3,_MAXATOMS_)
! phonon displacement distance
real(8) deltaph
! original lattice vectors
real(8) avec0(3,3)
! original inverse of lattice vector matrix
real(8) ainv0(3,3)
! original number of atoms
integer natoms0(_MAXSPECIES_)
integer natmtot0
! original atomic positions in Cartesian coordinates
real(8) atposc0(3,_MAXATOMS_,_MAXSPECIES_)
! original G-vector grid sizes
integer ngrid0(3)
integer ngrtot0
! original effective potentials
real(8), allocatable :: veffmt0(:,:,:)
real(8), allocatable :: veffir0(:)
! number of vectors for writing out frequencies and eigenvectors
integer nphwrt
! vectors in lattice coordinates for writing out frequencies and eigenvectors
real(8), allocatable :: vqlwrt(:,:)
! Coulomb pseudopotential
real(8) mustar
end module
