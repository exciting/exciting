!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: vecplot
! !INTERFACE:
!
!
Subroutine vecplot
! !DESCRIPTION:
      Use modinput
      use modplotlabels
      use modmpi, only : rank
!   Outputs a 2D or 3D vector field for plotting. The vector field can be the
!   magnetisation vector field, ${\bf m}$; the exchange-correlation magnetic
!   field, ${\bf B}_{\rm xc}$; or the electric field
!   ${\bf E}\equiv-\nabla V_{\rm C}$. The magnetisation is obtained from the
!   spin density matrix, $\rho_{\alpha\beta}$, by solving
!   $$ \rho_{\alpha\beta}({\bf r})=\frac{1}{2}\left(n({\bf r})
!    \delta_{\alpha\beta}+\sigma\cdot {\bf m}({\bf r})\right), $$
!   where $n\equiv\tr\rho_{\alpha\beta}$ is the total density. In the case of 2D
!   plots, the magnetisation vectors are still 3D, but are in the coordinate
!   system of the plane.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD)
!   Included electric field plots, August 2006 (JKD)
!EOP
!BOC
      Use modmain
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, lm
      Real (8) :: vl1 (3), vl2 (3), vc1 (3), vc2 (3), vc3 (3), vc4 (3), &
     & t1
     type(plotlabels),pointer ::labels
     type(plot2d_type),pointer::plot2ddef
     type(plot3d_type),pointer::plot3ddef
! allocatable arrays
      Real (8), Allocatable :: rvfmt (:, :, :, :)
      Real (8), Allocatable :: rvfir (:, :)
      If ((task .Eq. 72) .Or. (task .Eq. 73) .Or. &
      &   (task .Eq. 82) .Or. (task .Eq. 83)) Then
         If (.not.associated(input%groundstate%spin)) Then
            if (rank==0) then
              write (*,*)
              write (*, '("Error(vecplot): spin-unpolarised magnetisation field is zero")')
              Write (*,*)
              Stop
            end if
         End If
      End If
! initialise universal variables
      Call init0
! read magnetisation from file
      Call readstate
      Allocate (rvfmt(lmmaxvr, nrmtmax, natmtot, 3))
      Allocate (rvfir(ngrtot, 3))
      Select Case (task)
      Case (72, 73)
! magnetisation
         If (ncmag) Then
! non-collinear
            rvfmt (:, :, :, :) = magmt (:, :, :, :)
            rvfir (:, :) = magir (:, :)
         Else
! collinear
            rvfmt (:, :, :, 1:2) = 0.d0
            rvfir (:, 1:2) = 0.d0
            rvfmt (:, :, :, 3) = magmt (:, :, :, 1)
            rvfir (:, 3) = magir (:, 1)
         End If
      Case (82, 83)
! effective magnetic field
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
      Case (142, 143)
! electric field
         Call gradrf (vclmt, vclir, rvfmt, rvfir)
! use the negative of the gradient
         rvfmt (:, :, :, :) = - rvfmt (:, :, :, :)
         rvfir (:, :) = - rvfir (:, :)
      Case (152, 153)
         If (ndmag .Lt. 3) Then
            if (rank==0) then
              Write (*,*)
              Write (*, '("Error(vecplot): collinear m(r)xB(r) is zero")')
              Write (*,*)
              Stop
            end if
         End If
         Call rvfcross (magmt, bxcmt, magir, bxcir, rvfmt, rvfir)
      End Select
      Select Case (task)
      Case (72, 82, 142, 152)
! determine the projection of the magnetisation/field onto the plotting plane
         vl1 (:) = vclp2d (:, 2) - vclp2d (:, 1)
         vl2 (:) = vclp2d (:, 3) - vclp2d (:, 1)
         Call r3mv (input%structure%crystal%basevect, vl1, vc1)
         Call r3mv (input%structure%crystal%basevect, vl2, vc2)
         t1 = Sqrt (vc1(1)**2+vc1(2)**2+vc1(3)**2)
         vc1 (:) = vc1 (:) / t1
         t1 = Sqrt (vc2(1)**2+vc2(2)**2+vc2(3)**2)
         vc2 (:) = vc2 (:) / t1
         Call r3cross (vc1, vc2, vc3)
         t1 = Sqrt (vc3(1)**2+vc3(2)**2+vc3(3)**2)
         vc3 (:) = vc3 (:) / t1
         Call r3cross (vc3, vc1, vc2)
         t1 = Sqrt (vc2(1)**2+vc2(2)**2+vc2(3)**2)
         vc2 (:) = vc2 (:) / t1
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is)
                  Do lm = 1, lmmaxvr
                     vc4 (:) = rvfmt (lm, ir, ias, :)
                     rvfmt (lm, ir, ias, 1) = dot_product (vc4(:), &
                    & vc1(:))
                     rvfmt (lm, ir, ias, 2) = dot_product (vc4(:), &
                    & vc2(:))
                     rvfmt (lm, ir, ias, 3) = dot_product (vc4(:), &
                    & vc3(:))
                  End Do
               End Do
            End Do
         End Do
         Do ir = 1, ngrtot
            vc4 (:) = rvfir (ir, :)
            rvfir (ir, 1) = dot_product (vc4(:), vc1(:))
            rvfir (ir, 2) = dot_product (vc4(:), vc2(:))
            rvfir (ir, 3) = dot_product (vc4(:), vc3(:))
         End Do
         If (task .Eq. 72) Then
           plot2ddef=>input%properties%mvecfield%plot2d
           labels=>create_plotlablels("MAG","MAG2D",2)
		       call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"Magnetization","????","graceunit")
         Else If (task .Eq. 82) Then
           labels=>create_plotlablels("BXC","MAG2D",2)
           call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"BXC","????","graceunit")
		       plot2ddef=>input%properties%xcmvecfield%plot2d
         Else If (task .Eq. 142) Then
           labels=>create_plotlablels("EF","EF2D",2)
		       call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"EF2","????","graceunit")
           plot2ddef=>input%properties%electricfield%plot2d
         Else
           if (rank==0) write(*,*) "error in vecplot 2d selection"
         End If
         Call plot2d (labels, 3, input%groundstate%lmaxvr, lmmaxvr, rvfmt, &
         &            rvfir, plot2ddef)
         call destroy_plotlablels(labels)
         if (rank==0) then
           Write (*,*)
           Write (*, '("Info(vecplot):")')
           If (task .Eq. 72) Then
             Write (*, '(" 2D magnetisation density written to MAG2D.xml")')
           Else If (task .Eq. 82) Then
             Write (*, '(" 2D exchange-correlation field written to BXC2D.xml")')
           Else If (task .Eq. 142) Then
             Write (*, '(" 2D electric field written to EF2D.xml")')
           Else
             Write (*, '(" 2D m(r) x B_xc(r) written to MCBXC2D.xml")')
           End If
           Write (*,*)
         end if

      Case (73, 83, 143, 153)
         If (task .Eq. 73) Then
           labels=>create_plotlablels("MAG","MAG3D",3)
		       call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,4,"Magnetization","????","graceunit")
           plot3ddef=>input%properties%mvecfield%plot3d
         Else If (task .Eq. 83) Then
           labels=>create_plotlablels("BXC","BXC3D",3)
	         call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,4,"BXC","????","graceunit")
           plot3ddef=>input%properties%xcmvecfield%plot3d
         Else If (task .Eq. 143) Then
           labels=>create_plotlablels("EF","EF3D",3)
		       call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,4,"EF","????","graceunit")
		       plot3ddef=>input%properties%electricfield%plot3d
         Else
           labels=>create_plotlablels("MCBXC","MCBXC3D",3)
           call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
		       call set_plotlabel_axis(labels,4,"MCBXC","????","graceunit")
         End If
         Call plot3d (labels, 3, input%groundstate%lmaxvr, lmmaxvr, rvfmt, &
         &            rvfir,plot3ddef)
         call destroy_plotlablels(labels)
         if (rank==0) then
           Write (*,*)
           Write (*, '("Info(vecplot):")')
           If (task .Eq. 73) Then
             Write (*, '(" 3D magnetisation density written to MAG3D.xml")')
           Else If (task .Eq. 83) Then
             Write (*, '(" 3D exchange-correlation field written to BXC3D.xml")')
           Else If (task .Eq. 143) Then
             Write (*, '(" 3D electric field written to EF3D.xml")')
           Else
             Write (*, '(" 3D m(r) x B_xc(r) written to MCBXC3D.xml")')
           End If
           Write (*,*)
          end if
      End Select
      Deallocate (rvfmt, rvfir)
      Return
End Subroutine
!EOC
