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
      use modmpi, only : rank
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
        call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
        call set_plotlabel_axis(labels,2,"Potential","E_h/(ea_0)","graceunit")
        Call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & vclmt, vclir,input%properties%exccplot%plot1d)
        call destroy_plotlablels(labels)

  	    labels=>create_plotlablels("Potential","VXC1D",1)
	      call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
	      call set_plotlabel_axis(labels,2,"Exchange Correlation Potential","E_h/(ea_0)","graceunit")
        call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        &            vxcmt, vxcir, input%properties%exccplot%plot1d)
        call destroy_plotlablels(labels)
        if (rank==0) then
          Write (*,*)
          Write (*, '("Info(potplot):")')
          Write (*, '(" 1D Coulomb potential plot written to VCL1D.xml")')
          Write (*, '(" 1D exchange-correlation potential plot written to VXC1D.xml")')
          Write (*,*)
        end if
      endif ! 1d

      If (associated(input%properties%exccplot%plot2d)) then
        labels=>create_plotlablels("Potential","VCL2D",2)
	      call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,3,"Potential","E_h/(ea_0)","graceunit")
        call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vclmt, &
        &            vclir,input%properties%exccplot%plot2d)
        call destroy_plotlablels(labels)

        labels=>create_plotlablels("Potential","VXC2D",2)
        call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,3,"Potential","E_h/(ea_0)","graceunit")
        call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vxcmt, &
        &            vxcir,input%properties%exccplot%plot2d)
        call destroy_plotlablels(labels)
        if (rank==0) then
          Write (*,*)
          Write (*, '("Info(potplot):")')
          Write (*, '(" 2D Coulomb potential plot written to VCL2D.xml")')
          Write (*, '(" 2D exchange-correlation potential plot written to VXC2D.xml")')
          Write (*,*)
        end if
      endif ! 2d

      If (associated(input%properties%exccplot%plot3d)) then
        labels=>create_plotlablels("Potential","VCL3D",3)
	      call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,4,"Potential","E_h/(ea_0)","graceunit")
        call plot3d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vclmt, &
        &            vclir,input%properties%exccplot%plot3d)
        call destroy_plotlablels(labels)

        labels=>create_plotlablels("Potential","VXC3D",3)
        call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
        call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
        call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
	      call set_plotlabel_axis(labels,4,"Potential","E_h/(ea_0)","graceunit")
        call plot3d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, vxcmt, &
        &            vxcir,input%properties%exccplot%plot3d)
        call destroy_plotlablels(labels)
        if (rank==0) then
          Write (*,*)
          Write (*, '("Info(potplot):")')
          Write (*, '(" 3D Coulomb potential plot written to VCL3D.xml")')
          Write (*, '(" 3D exchange-correlation potential plot written to VXC3D.xml")')
          Write (*,*)
        end if
      End if ! 3d
      
      Return
End Subroutine
!EOC
