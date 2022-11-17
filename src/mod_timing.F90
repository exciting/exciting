
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!> Timing variables 
Module mod_timing
      implicit none 
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

      public :: stopwatch

contains

      !> Sirius timer.
      !>
      !> Note, if sirius_start_timer and sirius_stop_timer are straightforward,
      !> one could move them to exciting and general `stopwatch`, such that
      !> it can be used to time any code block.
      subroutine stopwatch(label, s)
#ifdef SIRIUS
      use sirius, only: sirius_start_timer, sirius_stop_timer
#endif
      !> Subroutine, module or code block label
      character(len=*), intent(in) :: label
      !> Start/stop integer: 0 or 1, respectively
      integer, intent(in) :: s
#ifdef SIRIUS
      if (s.eq.1) call sirius_start_timer(label)
      if (s.eq.0) call sirius_stop_timer(label)
#endif
      end subroutine stopwatch
End Module
