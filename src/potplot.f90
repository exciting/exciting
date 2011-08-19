!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: potplot
! !INTERFACE:
!
!
Subroutine potplot
! !USES:
      Use modinput
      Use modmain
      use modplotlabels
! !DESCRIPTION:
!   Outputs the exchange, correlation and Coulomb potentials, read in from
!   {\tt STATE.OUT}, for 1D, 2D or 3D plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! initialise universal variables
type(plotlabels),pointer ::labels
      Call init0
! read the density and potentials from file
      Call readstate
! write the potential plots to file
     If (associated(input%properties%exccplot%plot1d)) then

        labels=>create_plotlablels("Potential","VCL1D",1)
		 call set_plotlabel_axis(labels,1,"Distance","a_0")
		 call set_plotlabel_axis(labels,2,"Potential","E_h/(ea_0)")
         Call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & vclmt, vclir,input%properties%exccplot%plot1d)
         call destroy_plotlablels(labels)

  		  labels=>create_plotlablels("Potential","VXC1D",1)
		 call set_plotlabel_axis(labels,1,"Distance","a_0")
		 call set_plotlabel_axis(labels,2,"Exchange Correlation Potential","E_h/(ea_0)")
         Call plot1d ( labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & vxcmt, vxcir, input%properties%exccplot%plot1d)
         call destroy_plotlablels(labels)

         Write (*,*)
         Write (*, '("Info(potplot):")')
         Write (*, '(" 1D Coulomb potential plot written to VCL1D.xml")&
        &')
         Write (*, '(" 1D exchange-correlation potential plot written t&
        &o VXC1D.xml")')

    endif
     If (associated(input%properties%exccplot%plot2d)) then
          labels=>create_plotlablels("Potential","VCL2d",2)
		 call set_plotlabel_axis(labels,1,"a","lattice coordinate")
		 call set_plotlabel_axis(labels,2,"b","lattice coordinate")
		 call set_plotlabel_axis(labels,3,"Potential","E_h/(ea_0)")

         Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vclmt, &
        & vclir,input%properties%exccplot%plot2d)
        call destroy_plotlablels(labels)
         labels=>create_plotlablels("Potential","VXC2d",2)
		 call set_plotlabel_axis(labels,1,"a","lattice coordinate")
		 call set_plotlabel_axis(labels,2,"b","lattice coordinate")
		 call set_plotlabel_axis(labels,3,"Potential","E_h/(ea_0)")
         Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vxcmt, &
        & vxcir,input%properties%exccplot%plot2d)
         call destroy_plotlablels(labels)
         Write (*,*)
         Write (*, '("Info(potplot):")')
         Write (*, '(" 2D Coulomb potential plot written to VCL2d.xml")&
        &')
         Write (*, '(" 2D exchange-correlation potential plot written t&
        &o VXC2d.xml")')
     endif
      If (associated(input%properties%exccplot%plot3d)) then

         labels=>create_plotlablels("Potential","VCL3d",3)
		 call set_plotlabel_axis(labels,1,"a","lattice coordinate")
		 call set_plotlabel_axis(labels,2,"b","lattice coordinate")
		  call set_plotlabel_axis(labels,3,"c","lattice coordinate")
		 call set_plotlabel_axis(labels,4,"Potential","E_h/(ea_0)")
         Call plot3d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vclmt, &
        & vclir,input%properties%exccplot%plot3d)
         call destroy_plotlablels(labels)

           labels=>create_plotlablels("Potential","VXC3d",3)
		 call set_plotlabel_axis(labels,1,"a","lattice coordinate")
		 call set_plotlabel_axis(labels,2,"b","lattice coordinate")
		  call set_plotlabel_axis(labels,3,"c","lattice coordinate")
		 call set_plotlabel_axis(labels,4,"Potential","E_h/(ea_0)")
         Call plot3d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vxcmt, &
        & vxcir,input%properties%exccplot%plot3d)
         Close (50)
         Write (*,*)
         Write (*, '("Info(potplot):")')
         Write (*, '(" 3D Coulomb potential plot written to VCL3d.xml")&
        &')
         Write (*, '(" 3D exchange-correlation potential plot written t&
        &o VXC3d.xml")')
      End if
      Write (*,*)
      Return
End Subroutine
!EOC
