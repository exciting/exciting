
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
#include "maxdefinitions.inc"
Module mod_atoms
!
!--------------------------!
!     atomic variables     !
!--------------------------!
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
      Integer :: idxas (_MAXATOMS_, _MAXSPECIES_)
! molecule is .true. is the system is an isolated molecule
!replaced by inputstructure!replaced by inputstructurelogical::molecule
! primcell is .true. if primitive unit cell is to be found automatically
!replaced by inputstructure!replaced by inputstructurelogical::primcell
! atomic positions in lattice coordinates
!replaced by inputstructure!replaced by inputstructurereal(8)::atposl(3, _MAXATOMS_, _MAXSPECIES_)
! atomic positions in Cartesian coordinates
      Real (8) :: atposc (3, _MAXATOMS_, _MAXSPECIES_)
!
!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
!replaced by inputstructure!replaced by inputstructurecharacter(256)::sppath
! species filenames
!replaced by inputstructure!replaced by inputstructurecharacter(256)::spfname(_MAXSPECIES_)
! species name
      Character (256) :: spname (_MAXSPECIES_)
! species symbol
      Character (64) :: spsymb (_MAXSPECIES_)
! species nuclear charge
      Real (8) :: spzn (_MAXSPECIES_)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
!replaced by inputstructurelogical::ptnucl
! species electronic charge
      Real (8) :: spze (_MAXSPECIES_)
! species mass
      Real (8) :: spmass (_MAXSPECIES_)
! smallest radial point for each species
      Real (8) :: sprmin (_MAXSPECIES_)
! effective infinity for species
      Real (8) :: sprmax (_MAXSPECIES_)
! number of radial points to effective infinity for each species
      Integer :: spnr (_MAXSPECIES_)
! maximum spnr over all the species
      Integer :: spnrmax
! maximum allowed states for each species
      Integer, Parameter :: maxspst = 40
! number of states for each species
      Integer :: spnst (_MAXSPECIES_)
! maximum spnst over all the species
      Integer :: spnstmax
! state principle quantum number for each species
      Integer :: spn (maxspst, _MAXSPECIES_)
! state l value for each species
      Integer :: spl (maxspst, _MAXSPECIES_)
! state k value for each species
      Integer :: spk (maxspst, _MAXSPECIES_)
! spcore is .true. if species state is core
      Logical :: spcore (maxspst, _MAXSPECIES_)
! state eigenvalue for each species
      Real (8) :: speval (maxspst, _MAXSPECIES_)
! state occupancy for each species
      Real (8) :: spocc (maxspst, _MAXSPECIES_)
! species radial mesh
      Real (8), Allocatable :: spr (:, :)
! species charge density
      Real (8), Allocatable :: sprho (:, :)
! species self-consistent potential
      Real (8), Allocatable :: spvr (:, :)
!
End Module
!
