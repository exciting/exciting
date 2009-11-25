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
