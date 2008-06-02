
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

integer function iplocnr(ivp,ngridp)
  ! locate p-point with index on grid following the convention that the
  ! first coordinate runs fastest.
  implicit none
  ! arguments
  integer, intent(in) :: ivp(3),ngridp(3)
  ! this should be consistent with ipmap of "genppts.f90"
  iplocnr = 1 + ivp(1) + ngridp(1)*ivp(2) + ngridp(1)*ngridp(2)*ivp(3)
end function iplocnr
