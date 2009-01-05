#include "maxdefinitions.inc"
module mod_force
!-------------------------!
!     force variables     !
!-------------------------!
! tforce is .true. if force should be calculated
logical tforce
! tfibs is .true. if the IBS contribution to the force is to be calculated
logical tfibs
! Hellmann-Feynman force on each atom
real(8), allocatable :: forcehf(:,:)
! core correction to force on each atom
real(8), allocatable :: forcecr(:,:)
! IBS core force on each atom
real(8), allocatable :: forceibs(:,:)
! total force on each atom
real(8), allocatable :: forcetot(:,:)
! previous total force on each atom
real(8), allocatable :: forcetp(:,:)
! maximum force magnitude over all atoms
real(8) forcemax
! default step size parameter for structural optimisation
real(8) tau0atm
! step size parameters for each atom
real(8), allocatable :: tauatm(:)
end module
