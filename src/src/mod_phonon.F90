
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_phonon
!--------------------------!
!     phonon variables     !
!--------------------------!
! number of primitive unit cells in phonon supercell
      Integer :: nphcell
! Cartesian offset vectors for each primitive cell in the supercell
      Real (8) :: vphcell (3, _MAXATOMS_)
! phonon displacement distance
!replaced by inputstructurereal(8)::deltaph
! original lattice vectors
      Real (8) :: avec0 (3, 3)
! original inverse of lattice vector matrix
      Real (8) :: ainv0 (3, 3)
! original number of atoms
      Integer :: natoms0 (_MAXSPECIES_)
      Integer :: natmtot0
! original atomic positions in Cartesian coordinates
      Real (8) :: atposc0 (3, _MAXATOMS_, _MAXSPECIES_)
! original G-vector grid sizes
      Integer :: ngrid0 (3)
      Integer :: ngrtot0
! original effective potentials
      Real (8), Allocatable :: veffmt0 (:, :, :)
      Real (8), Allocatable :: veffir0 (:)
! number of vectors for writing out frequencies and eigenvectors
      Integer :: nphwrt
! vectors in lattice coordinates for writing out frequencies and eigenvectors
      Real (8), Allocatable :: vqlwrt (:, :)
! Coulomb pseudopotential
!replaced by inputstructurereal(8)::mustar
! file suffix for DYN files
      Character(256) :: filextdyn
! subdirectory name for phonon calculation
      Character(256) :: phdirname
End Module
!
