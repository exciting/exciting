
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_corestate
!
!------------------------------!
!     core state variables     !
!------------------------------!
! eigenvalues for core states
      Real (8), Allocatable :: evalcr (:, :)
! radial wavefunctions for core states
      Real (8), Allocatable :: rwfcr (:, :, :, :)
! radial charge density for core states
      Real (8), Allocatable :: rhocr (:, :)
! true for a frozen-core calculation (core state wavefunctions, densities and
! energies calculated only in the first iteration)
!replaced by inputstructurelogical :: frozencore
End Module
!
