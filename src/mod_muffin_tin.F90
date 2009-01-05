#include "maxdefinitions.inc"
module mod_muffin_tin
!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! radial function integration and differentiation polynomial order
integer nprad
! number of muffin-tin radial points for each species
integer nrmt(_MAXSPECIES_)
! maximum nrmt over all the species
integer nrmtmax
! autormt is .true. for automatic determination of muffin-tin radii
logical autormt
! parameters for determining muffin-tin radii automatically
real(8) rmtapm(2)
! muffin-tin radii
real(8) rmt(_MAXSPECIES_)
! species for which the muffin-tin radius will be used for calculating gkmax
integer isgkmax
! radial step length for coarse mesh
integer lradstp
! number of coarse radial mesh points
integer nrcmt(_MAXSPECIES_)
! maximum nrcmt over all the species
integer nrcmtmax
! coarse muffin-tin radial mesh
real(8), allocatable :: rcmt(:,:)
! maximum allowable angular momentum for augmented plane waves

! maximum angular momentum for augmented plane waves
integer lmaxapw
! (lmaxapw+1)^2
integer lmmaxapw
! maximum angular momentum for potentials and densities
integer lmaxvr
! (lmaxvr+1)^2
integer lmmaxvr
! maximum angular momentum used when evaluating the Hamiltonian matrix elements
integer lmaxmat
! (lmaxmat+1)^2
integer lmmaxmat
! fraction of muffin-tin radius which constitutes the inner part
real(8) fracinr
! maximum angular momentum in the inner part of the muffin-int
integer lmaxinr
! (lmaxinr+1)^2
integer lmmaxinr
! number of radial points to the inner part of the muffin-tin
integer nrmtinr(_MAXSPECIES_)
! index to (l,m) pairs
integer, allocatable :: idxlm(:,:)

end module
