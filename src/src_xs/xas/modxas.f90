! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE: modxas
! !DESCRIPTION:
!   Contains additional global variables required for XAS calculations within the BSE EXCITING code.
!
! !REVISION HISTORY:
!
!   Created JUNE 2015 by Christian Vorwerk
!EOP   
!BOC
#include "maxdefinitions.inc"
module modxas
!----------------------------!
!     Core states            !
!----------------------------!
! WARNING: the definitions here are not identical to the one in modgw!
! See xasinit.f90 for how core treatment is done in XAS-calculations
! maximum number of core states per atom
      integer(4) :: ncmax
! Max. num of core states including lm
      integer(4) :: nclm 
! Max. L of core states
      integer(4) :: lcoremax      
! Total num of core states over all atoms including lm
      integer(4) :: ncg  
! number of core states of given species and atom
      integer(4) :: ncore  
! number of initial states for XAS calculation
       integer(4) :: nxas
! radial wavefunctions for given kappa and m
      real(8), allocatable :: ucore(:,:)
!  core energies for given kappa and m
	  real(8), allocatable :: ecore(:)
!  corresponding ml for a given (l,j,mj)-state
	  integer(4), allocatable :: mj2ml(:,:)
! prefactor for the spherical harmonics within the spin spherical harmonics
	  real(8), allocatable :: preml(:,:)
! unique lowest and highest core state for XAS calculation
	  integer(4) :: xasstart, xasstop
! l-value of the considered XAS edge
	  integer(4) :: lxas
!----------------------------!
!     general                !
!----------------------------!
! shortcut for atomic position array
      real(8) :: atposl(3,_MAXATOMS_,_MAXSPECIES_)
! muffin-tin volume, relative to cell volume
      real(8) :: vmt(_MAXSPECIES_)
! plane wave matrix elements array (o-o part)
      Complex (8), Allocatable :: xioo (:, :, :)
! plane wave matrix elements array (u-u part)
      Complex (8), Allocatable :: xiuu (:, :, :)

end module modxas
!EOC
