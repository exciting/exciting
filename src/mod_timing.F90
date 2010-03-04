
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

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
