
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: iplocnr
! !INTERFACE:
integer function iplocnr(ivp,ngridp)

use modmain, only: iqmap
use modxs

! !INPUT/OUTPUT PARAMETERS:
!   ivp    : integer coordinates of p-point (in,integer(3))
!   ngridp : p-point set grid sizes (in,integer(3))
! !DESCRIPTION:
!   Locate {\bf p}-point with index on grid following the convention that the
!   first coordinate runs fastest.
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: ivp(3),ngridp(3)
  ! this should be consistent with ipmap of "genppts.f90"
  iplocnr = 1 + ivp(1) + ngridp(1)*ivp(2) + ngridp(1)*ngridp(2)*ivp(3)

!!!iplocnr=iqmap(ivp(1),ivp(2),ivp(3))

end function iplocnr
!EOC
