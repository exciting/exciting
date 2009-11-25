!
!
#include "maxdefinitions.inc"
Module mod_timing
!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
      Real (8) :: timeinit
! Hamiltonian and overlap matrix set up
      Real (8) :: timemat
! first-variational calculation
      Real (8) :: timefv
! second-variational calculation
      Real (8) :: timesv
! charge density calculation
      Real (8) :: timerho
! potential calculation
      Real (8) :: timepot
! force calculation
      Real (8) :: timefor
End Module
!
