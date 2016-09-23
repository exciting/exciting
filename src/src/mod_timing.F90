
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
! input and output
      Real (8) :: timeio
! muffin-tin operations
      Real (8) :: timemt
! mixer
      Real (8) :: timemixer
! matching coefficients 
      Real (8) :: timematch

! Parts of Hamiltonian and overlap matrix set up
      Real (8) :: time_hmlaan
      Real (8) :: time_hmlalon
      Real (8) :: time_hmllolon
      Real (8) :: time_olpaan
      Real (8) :: time_olpalon
      Real (8) :: time_olplolon
      Real (8) :: time_hmlistln
      Real (8) :: time_olpistln

! Parts of initialisation
      Real (8) :: time_init0
      Real (8) :: time_init1
      Real (8) :: time_density_init
      Real (8) :: time_pot_init

! Radial solvers   
      Real (8) :: time_rdirac
      Real (8) :: time_rschrod 
! OEP 
      Real (8) :: time_oepvnl
      Real (8) :: time_oep_iter

End Module
!
