!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhoplot
! !INTERFACE:
!
!
Subroutine rhoplot
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Outputs the charge density, read in from {\tt STATE.OUT}, for 1D, 2D or 3D
!   plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! initialise universal variables
      Call init0
! read density from file
      Call readstate
! write the density plot to file
      If (associated(input%properties%chargedensityplot%plot1d)) Then
!
         Call plot1d ("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot1d)
!
         Write (*,*)
         Write (*, '("Info(rhoplot):")')
         Write (*, '(" 1D density plot written to RHO1D.OUT")')
         Write (*, '(" vertex location lines written to RHOLINES.OUT")')
      End If
      If (associated(input%properties%chargedensityplot%plot2d)) Then
!
         Call plot2d ("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot2d)
!
         Write (*,*)
         Write (*, '("Info(rhoplot): 2D density plot written to RHO2D.O&
        &UT")')
      End If
      If (associated(input%properties%chargedensityplot%plot3d)) Then
         Call plot3d ("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot3d)
         Write (*,*)
         Write (*, '("Info(rhoplot): 3D density plot written to RHO3D.O&
        &UT")')
      End If
      Write (*,*)
      Return
End Subroutine
!EOC
