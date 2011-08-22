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
      use modplotlabels

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
      type(plotlabels),pointer ::labels

      Call init0
! read density from file
      Call readstate
! write the density plot to file
      If (associated(input%properties%chargedensityplot%plot1d)) Then
        labels=>create_plotlablels("RHO","RHO",1)
		 call set_plotlabel_axis(labels,1,"Distance","a_0")
		 call set_plotlabel_axis(labels,2,"Density","???")
         Call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot1d)
         call destroy_plotlablels(labels)
         Write (*,*)
         Write (*, '("Info(rhoplot):")')
         Write (*, '(" 1D density plot written to RHO1D.xml")')

      End If
      If (associated(input%properties%chargedensityplot%plot2d)) Then
!
		 labels=>create_plotlablels("RHO","RHO",2)
		 call set_plotlabel_axis(labels,1,"a","1")
		 call set_plotlabel_axis(labels,2,"b","1")
		 call set_plotlabel_axis(labels,3,"Density","???")
         Call plot2d ( labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot2d)
         call destroy_plotlablels(labels)
         Write (*,*)
         Write (*, '("Info(rhoplot): 2D density plot written to RHO2D.xml")')
      End If
      If (associated(input%properties%chargedensityplot%plot3d)) Then
         labels=>create_plotlablels("RHO","RHO",3)
		 call set_plotlabel_axis(labels,1,"a","1")
		 call set_plotlabel_axis(labels,2,"b","1")
		 call set_plotlabel_axis(labels,3,"b","1")
		 call set_plotlabel_axis(labels,4,"Density","???")
         Call plot3d ( labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot3d)
         call destroy_plotlablels(labels)
         Write (*,*)
         Write (*, '("Info(rhoplot): 3D density plot written to RHO3D.O&
        &UT")')
      End If
      Write (*,*)
      Return
End Subroutine
!EOC
