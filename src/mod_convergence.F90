
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_convergence
!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
!replaced by inputstructureinteger::maxscl
! current self-consistent loop number
      Integer :: iscl
! current structure optimization step
      Integer :: istep
! current structure optimization step
      Logical :: lstep
! effective potential convergence tolerance
!replaced by inputstructurereal(8)::epspot
! energy convergence tolerance
!replaced by inputstructurereal(8)::epsengy
! force convergence tolerance
!replaced by inputstructurereal(8)::epsforce
!curent convergence
      Real (8) :: currentconvergence
End Module
!
