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
end module
