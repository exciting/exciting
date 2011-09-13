!
!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dbxcplot
      Use modmain
      Use modinput
      use modplotlabels
      Implicit None
! local variables
      Integer :: idm, is, ia, ias, ir
      type(plotlabels),pointer ::labels
! allocatable arrays
      Real (8), Allocatable :: rvfmt (:, :, :, :)
      Real (8), Allocatable :: rvfir (:, :)
      Real (8), Allocatable :: rfmt (:, :, :)
      Real (8), Allocatable :: rfir (:)
      Real (8), Allocatable :: grfmt (:, :, :, :)
      Real (8), Allocatable :: grfir (:, :)
      If ( .Not. associated(input%groundstate%spin)) Then
         Write (*,*)
         Write (*, '("Error(dbxcplot): spin-unpolarised field is zero")&
        &')
         Write (*,*)
         Stop
      End If
! initialise universal variables
      Call init0
! read magnetisation from file
      Call readstate
      Allocate (rvfmt(lmmaxvr, nrmtmax, natmtot, 3))
      Allocate (rvfir(ngrtot, 3))
      Allocate (rfmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (rfir(ngrtot))
      Allocate (grfmt(lmmaxvr, nrmtmax, natmtot, 3))
      Allocate (grfir(ngrtot, 3))
      If (ncmag) Then
! non-collinear
         rvfmt (:, :, :, :) = bxcmt (:, :, :, :)
         rvfir (:, :) = bxcir (:, :)
      Else
! collinear
         rvfmt (:, :, :, 1:2) = 0.d0
         rvfir (:, 1:2) = 0.d0
         rvfmt (:, :, :, 3) = bxcmt (:, :, :, 1)
         rvfir (:, 3) = bxcir (:, 1)
      End If
      rfmt (:, :, :) = 0.d0
      rfir (:) = 0.d0
      Do idm = 1, 3
         Call gradrf (rvfmt(:, :, :, idm), rvfir(:, idm), grfmt, grfir)
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is)
                  rfmt (:, ir, ias) = rfmt (:, ir, ias) + grfmt (:, ir, &
                 & ias, idm)
               End Do
            End Do
         End Do
         rfir (:) = rfir (:) + grfir (:, idm)
      End Do
      If (associated(input%properties%gradmvecfield%plot1d)) Then
         labels=>create_plotlablels("DBXC","DBXC",1)
		 call set_plotlabel_axis(labels,1,"Distance","a_0","graceunit")
		 call set_plotlabel_axis(labels,2,"gradmvecfield","???","graceunit")
         Call plot1d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rfmt, rfir, input%properties%gradmvecfield%plot1d)
         call destroy_plotlablels(labels)
         Write (*,*)
         Write (*, '("Info(dbxcplot):")')
         Write (*, '(" 1D divergence of exchange-correlation field writ&
        &ten to DBXC1D.xml")')


      End If
      If (associated(input%properties%gradmvecfield%plot2d)) Then
         labels=>create_plotlablels("DBXC","DBXC",2)
		 call set_plotlabel_axis(labels,1,"a","1","graceunit")
		 call set_plotlabel_axis(labels,2,"b","1","graceunit")
		 call set_plotlabel_axis(labels,3,"gradmvecfield","???","graceunit")
         Call plot2d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rfmt, rfir, input%properties%gradmvecfield%plot2d)
         call destroy_plotlablels(labels)
         Write (*, '("Info(dbxcplot):")')
         Write (*, '(" 2D divergence of exchange-correlation field writ&
        &ten to DBXC2d.xml")')
      End If
      If (associated(input%properties%gradmvecfield%plot3d)) Then
         labels=>create_plotlablels("DBXC","DBXC",3)
		 call set_plotlabel_axis(labels,1,"a","1","graceunit")
		 call set_plotlabel_axis(labels,2,"b","1","graceunit")
		 call set_plotlabel_axis(labels,3,"b","1","graceunit")
		 call set_plotlabel_axis(labels,4,"gradmvecfield","???","graceunit")
         Call plot3d (labels, 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rfmt, rfir, input%properties%gradmvecfield%plot3d)
        call destroy_plotlablels(labels)
         Write (*, '("Info(dbxcplot):")')
         Write (*, '(" 3D divergence of exchange-correlation field writ&
        &ten to DBXC3d.xml")')
      End If
      Write (*,*)
      Deallocate (rvfmt, rvfir, rfmt, rfir, grfmt, grfir)
      Return
End Subroutine
