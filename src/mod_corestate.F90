#include "maxdefinitions.inc"
module mod_corestate

!------------------------------!
!     core state variables     !
!------------------------------!
! eigenvalues for core states
real(8), allocatable :: evalcr(:,:)
! radial wavefunctions for core states
real(8), allocatable :: rwfcr(:,:,:,:)
! radial charge density for core states
real(8), allocatable :: rhocr(:,:)
! true for a frozen-core calculation (core state wavefunctions, densities and
! energies calculated only in the first iteration)
logical :: frozencore
end module
