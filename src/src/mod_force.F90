
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_force
!-------------------------!
!     force variables     !
!-------------------------!
! tforce is .true. if force should be calculated
!replaced by inputstructurelogical::tforce
! tfibs is .true. if the IBS contribution to the force is to be calculated
!replaced by inputstructurelogical::tfibs
! Hellmann-Feynman force on each atom
      Real (8), Allocatable :: forcehf (:, :)
! core correction to force on each atom
      Real (8), Allocatable :: forcecr (:, :)
! IBS core force on each atom
      Real (8), Allocatable :: forceibs (:, :)
! dispersion correction force on each atom
      Real (8), Allocatable :: force_disp (:, :)
! total force on each atom
      Real (8), Allocatable :: forcetot (:, :)
! previous total force on each atom
      Real (8), Allocatable :: forcetp (:, :)
! maximum force magnitude over all atoms
      Real (8) :: forcemax
! default step size parameter for structural optimisation
!replaced by inputstructurereal(8)::taunewton
! step size parameters for each atom
      Real (8), Allocatable :: tauatm (:)
      Real (8), Allocatable :: tauxyz (:, :)
End Module
!
