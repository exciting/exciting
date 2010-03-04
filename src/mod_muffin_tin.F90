
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_muffin_tin
!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! radial function integration and differentiation polynomial order
!replaced by inputstructureinteger::nprad
! number of muffin-tin radial points for each species
      Integer :: nrmt (_MAXSPECIES_)
! maximum nrmt over all the species
      Integer :: nrmtmax
! autormt is .true. for automatic determination of muffin-tin radii
!replaced by inputstructurelogical::autormt
! parameters for determining muffin-tin radii automatically
!replaced by inputstructurereal(8)::rmtapm(2)
! muffin-tin radii
      Real (8) :: rmt (_MAXSPECIES_)
! species for which the muffin-tin radius will be used for calculating gkmax
!replaced by inputstructureinteger::isgkmax
! radial step length for coarse mesh
!replaced by inputstructureinteger::lradstp
! number of coarse radial mesh points
      Integer :: nrcmt (_MAXSPECIES_)
! maximum nrcmt over all the species
      Integer :: nrcmtmax
! coarse muffin-tin radial mesh
      Real (8), Allocatable :: rcmt (:, :)
! maximum allowable angular momentum for augmented plane waves
!
! maximum angular momentum for augmented plane waves
!replaced by inputstructureinteger::lmaxapw
! (lmaxapw+1)^2
      Integer :: lmmaxapw
! maximum angular momentum for potentials and densities
!replaced by inputstructureinteger::lmaxvr
! (lmaxvr+1)^2
      Integer :: lmmaxvr
! maximum angular momentum used when evaluating the Hamiltonian matrix elements
!replaced by inputstructureinteger::lmaxmat
! (lmaxmat+1)^2
      Integer :: lmmaxmat
! fraction of muffin-tin radius which constitutes the inner part
!replaced by inputstructurereal(8)::fracinr
! maximum angular momentum in the inner part of the muffin-int
!replaced by inputstructureinteger::lmaxinr
! (lmaxinr+1)^2
      Integer :: lmmaxinr
! number of radial points to the inner part of the muffin-tin
      Integer :: nrmtinr (_MAXSPECIES_)
! index to (l,m) pairs
      Integer, Allocatable :: idxlm (:, :)
!------------------------------!
!     tolerance parameters     !
!------------------------------!
! energy convergence tolerance for Dirac equation solver
      Real (8) :: epsedirac
! potential convergence tolerance for atomistic calculation
      Real (8) :: epspotatom
End Module
!
