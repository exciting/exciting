
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE: modmain
! !DESCRIPTION:
!   Contains all the global variables required by the spacegroup code.
!
! !REVISION HISTORY:
!   Created October 2006 (JKD)
!EOP
!BOC
module modmain

!-------------------------------!
!     space group variables     !
!-------------------------------!
! Hermann-Mauguin symbol
character(20) hrmg
! space-group number
character(20) num
! Schoenflies symbol
character(20) schn
! Hall symbol
character(20) hall

!----------------------------!
!     lattice parameters     !
!----------------------------!
! number of unit cells
integer ncell(3)
! lattice vector lengths
real(8) a,b,c
! lattice vector angles
real(8) ab,ac,bc
! lattice vectors stored column-wise
real(8) avec(3,3)
! inverse of lattice vector matrix
real(8) ainv(3,3)
! any vector with length less than epslat is considered zero
real(8), parameter :: epslat=1.d-6

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=1000
! number of species
integer nspecies
! number of atoms for each species
integer natoms(maxspecies)
! total number of atoms
integer natmtot
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! maximum allowed Wyckoff positions
integer, parameter :: maxwpos=100
! number of Wyckoff positions
integer nwpos(maxspecies)
! Wyckoff positions
real(8) wpos(3,maxwpos,maxspecies)
! atomic positions in lattice coordinates
real(8) atposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc(3,maxatoms,maxspecies)
! magnetic fields
real(8) bfcmt(3,maxatoms,maxspecies)

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species file names
character(256) spfname(maxspecies)
! species name
character(256) spname(maxspecies)
! species symbol
character(256) spsymb(maxspecies)

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer version(3)
data version / 1,1,3 /

end module
!EOC
