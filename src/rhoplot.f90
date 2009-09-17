



! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoplot
! !INTERFACE:


subroutine rhoplot
! !USES:
use modinput
use modmain
! !DESCRIPTION:
!   Outputs the charge density, read in from {\tt STATE.OUT}, for 1D, 2D or 3D
!   plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! initialise universal variables
call init0
! read density from file
call readstate
! write the density plot to file
if(associated(input%properties%chargedensityplot%plot1d)) then

  call plot1d("RHO" ,1, input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir,input%properties%chargedensityplot%plot1d)

  write(*, *)
  write(*, '("Info(rhoplot):")')
  write(*, '(" 1D density plot written to RHO1D.OUT")')
  write(*, '(" vertex location lines written to RHOLINES.OUT")')
endif
if(associated(input%properties%chargedensityplot%plot2d)) then

  call plot2d("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir,input%properties%chargedensityplot%plot2d)

  write(*, *)
  write(*, '("Info(rhoplot): 2D density plot written to RHO2D.OUT")')
endif
if(associated(input%properties%chargedensityplot%plot3d)) then
  call plot3d("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, rhomt, rhoir,input%properties%chargedensityplot%plot3d)
  write(*, *)
  write(*, '("Info(rhoplot): 3D density plot written to RHO3D.OUT")')
endif
write(*, *)
return
end subroutine
!EOC
