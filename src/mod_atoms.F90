
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> atomic variables 
Module mod_atoms
      use constants, only: maxatoms, maxspecies
      implicit none 

! maximum allowed species
!
! maximum allowed atoms per species
!
! number of species
      Integer :: nspecies
! number of atoms for each species
      Integer, Allocatable :: natoms (:)
! maximum number of atoms over all the species
      Integer :: natmmax
! total number of atoms
      Integer :: natmtot
! index to atoms and species
      Integer :: idxas (maxatoms, maxspecies)
! molecule is .true. is the system is an isolated molecule
!replaced by inputstructure!replaced by inputstructurelogical::molecule
! primcell is .true. if primitive unit cell is to be found automatically
!replaced by inputstructure!replaced by inputstructurelogical::primcell
! atomic positions in lattice coordinates
!replaced by inputstructure!replaced by inputstructurereal(8)::atposl(3, maxatoms, maxspecies)
! atomic positions in Cartesian coordinates
      Real (8) :: atposc (3, maxatoms, maxspecies)
      Real (8) :: atsave (3, maxatoms, maxspecies)
      Real (8) :: atposcm (3, maxatoms, maxspecies)
      Real (8) :: atposc_1 (3, maxatoms, maxspecies)
      Real (8) :: tsec1, tsec2
!
!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
!replaced by inputstructure!replaced by inputstructurecharacter(256)::sppath
! species filenames
!replaced by inputstructure!replaced by inputstructurecharacter(256)::spfname(maxspecies)
! species name
      Character (256) :: spname (maxspecies)
! species symbol
      Character (64) :: spsymb (maxspecies)
! species nuclear charge
      Real (8) :: spzn (maxspecies)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
!replaced by inputstructurelogical::ptnucl
! species electronic charge
      Real (8) :: spze (maxspecies)
! species mass
      Real (8) :: spmass (maxspecies)
! smallest radial point for each species
      Real (8) :: sprmin (maxspecies)
! effective infinity for species
      Real (8) :: sprmax (maxspecies)
! number of radial points to effective infinity for each species
      Integer :: spnr (maxspecies)
! maximum spnr over all the species
      Integer :: spnrmax
! maximum allowed states for each species
      Integer, Parameter :: maxspst = 40
! number of states for each species
      Integer :: spnst (maxspecies)
! maximum spnst over all the species
      Integer :: spnstmax
! state principle quantum number for each species
      Integer :: spn (maxspst, maxspecies)
! state l value for each species
      Integer :: spl (maxspst, maxspecies)
! state k value for each species
      Integer :: spk (maxspst, maxspecies)
! spcore is .true. if species state is core
      Logical :: spcore (maxspst, maxspecies)
! state eigenvalue for each species
      Real (8) :: speval (maxspst, maxspecies)
! state occupancy for each species
      Real (8) :: spocc (maxspst, maxspecies)
! species radial mesh
      Real (8), Allocatable :: spr (:, :)
! species charge density
      Real (8), Allocatable :: sprho (:, :)
! species self-consistent potential
      Real (8), Allocatable :: spvr (:, :)
!
End Module
!
