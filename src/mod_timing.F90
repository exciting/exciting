

#include "maxdefinitions.inc"
module mod_timing
!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
real(8)::timeinit
! Hamiltonian and overlap matrix set up
real(8)::timemat
! first-variational calculation
real(8)::timefv
! second-variational calculation
real(8)::timesv
! charge density calculation
real(8)::timerho
! potential calculation
real(8)::timepot
! force calculation
real(8)::timefor
end module

