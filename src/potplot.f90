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
      Call init0
! read the density and potentials from file
      Call readstate
! write the potential plots to file
      Select Case (task)
      Case (41)
         Open (50, File='VCL1D.OUT', Action='WRITE', Form='FORMATTED')
         Open (51, File='VLINES.OUT', Action='WRITE', Form='FORMATTED')
         Call plot1d (50, 51, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & vclmt, vclir)
         Close (50)
         Close (51)
         Open (50, File='VXC1D.OUT', Action='WRITE', Form='FORMATTED')
         Open (51, File='VLINES.OUT', Action='WRITE', Form='FORMATTED')
         Call plot1d (50, 51, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & vxcmt, vxcir)
         Close (50)
         Close (51)
         Write (*,*)
         Write (*, '("Info(potplot):")')
         Write (*, '(" 1D Coulomb potential plot written to VCL1D.OUT")&
        &')
         Write (*, '(" 1D exchange-correlation potential plot written t&
        &o VXC1D.OUT")')
         Write (*, '(" vertex location lines written to VLINES.OUT")')
      Case (42)
         Open (50, File='VCL2d.xml', Action='WRITE', Form='FORMATTED')
         Call plot2d (50, 1, input%groundstate%lmaxvr, lmmaxvr, vclmt, &
        & vclir)
         Close (50)
         Open (50, File='VXC2d.xml', Action='WRITE', Form='FORMATTED')
         Call plot2d (50, 1, input%groundstate%lmaxvr, lmmaxvr, vxcmt, &
        & vxcir)
         Close (50)
         Write (*,*)
         Write (*, '("Info(potplot):")')
         Write (*, '(" 2D Coulomb potential plot written to VCL2d.xml")&
        &')
         Write (*, '(" 2D exchange-correlation potential plot written t&
        &o VXC2d.xml")')
      Case (43)
         Open (50, File='VCL3d.xml', Action='WRITE', Form='FORMATTED')
         Call plot3d (50, 1, input%groundstate%lmaxvr, lmmaxvr, vclmt, &
        & vclir)
         Close (50)
         Open (50, File='VXC3d.xml', Action='WRITE', Form='FORMATTED')
         Call plot3d (50, 1, input%groundstate%lmaxvr, lmmaxvr, vxcmt, &
        & vxcir)
         Close (50)
         Write (*,*)
         Write (*, '("Info(potplot):")')
         Write (*, '(" 3D Coulomb potential plot written to VCL3d.xml")&
        &')
         Write (*, '(" 3D exchange-correlation potential plot written t&
        &o VXC3d.xml")')
      End Select
      Write (*,*)
      Return
End Subroutine
!EOC
