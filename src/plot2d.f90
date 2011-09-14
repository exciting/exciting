! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: plot2d
! !INTERFACE:
!
!
Subroutine plot2d (labels, nf, lmax, ld, rfmt, rfir, plotdef)
! !USES:
      Use modinput
 use mod_muffin_tin
  use mod_atoms
  use mod_Gvector
      Use FoX_wxml
      use modmpi
  use modplotlabels
  use mod_plotting
! !INPUT/OUTPUT PARAMETERS:
!   fname : plot file name character(len=*)
!   nf   : number of functions (in,integer)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngrtot,nf))
!   plotdef:type(plot2d) defines plot region
! !DESCRIPTION:
!   Produces a 2D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} on the parallelogram defined by the corner vertices in the global
!   array {\tt vclp2d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      type(plotlabels), Intent (In) :: labels
      Integer, Intent (In) :: nf
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: rfmt (ld, nrmtmax, natmtot, nf)
      Real (8), Intent (In) :: rfir (ngrtot, nf)
      Type (plot2d_type) :: plotdef

! local variables
      Integer :: i, ip, ip1, ip2,ifunction
      Real (8) :: vl1 (3), vl2 (3), vc1 (3), vc2 (3),delta(3)
      Real (8) :: d1, d2, d12, t1, t2, t3, t4
      Character (128) :: buffer, buffer1
      character (20)::buffer20
      Type (xmlf_t), Save :: xf
! allocatable arrays
      Real (8), Allocatable :: vpl (:, :)
      Real (8), Allocatable :: fp (:, :)
!external functions
      Real(8),external::DNRM2
 If (rank .Eq. 0) Then
      write(buffer,*) labels%filename , "2D.XML"
      Call xml_OpenFile (trim(buffer), xf, replace=.True.,  pretty_print=.True.)
      Call xml_NewElement (xf, "plot2d")
!
      If ((nf .Lt. 1) .Or. (nf .Gt. 4)) Then
         Write (*,*)
         Write (*, '("Error(plot2d): invalid number of functions : ", I&
        &8)') nf
         Write (*,*)
         Stop
      End If
! allocate local arrays
      Allocate (vpl(3, &
     & plotdef%parallelogram%grid(1)*plotdef%parallelogram%grid(2)))
      Allocate &
     & (fp(plotdef%parallelogram%grid(1)*plotdef%parallelogram%grid(2), &
     & nf))
! generate 2D grid
      vl1 (:) = plotdef%parallelogram%pointarray(1)%point%coord - &
     & plotdef%parallelogram%origin%coord
      vl2 (:) = plotdef%parallelogram%pointarray(2)%point%coord - &
     & plotdef%parallelogram%origin%coord
      vc1 (:) = vl1 (1) * input%structure%crystal%basevect(:, 1) + vl1 &
     & (2) * input%structure%crystal%basevect(:, 2) + vl1 (3) * &
     & input%structure%crystal%basevect(:, 3)
      vc2 (:) = vl2 (1) * input%structure%crystal%basevect(:, 1) + vl2 &
     & (2) * input%structure%crystal%basevect(:, 2) + vl2 (3) * &
     & input%structure%crystal%basevect(:, 3)
      d1 = Sqrt (vc1(1)**2+vc1(2)**2+vc1(3)**2)
      d2 = Sqrt (vc2(1)**2+vc2(2)**2+vc2(3)**2)
      If ((d1 .Lt. input%structure%epslat) .Or. (d2 .Lt. &
     & input%structure%epslat)) Then
         Write (*,*)
         Write (*, '("Error(plot2d): zero length plotting vectors")')
         Write (*,*)
         Stop
      End If
      d12 = (vc1(1)*vc2(1)+vc1(2)*vc2(2)+vc1(3)*vc2(3)) / (d1*d2)
      ip = 0
      Do ip1 = 0, plotdef%parallelogram%grid(1) - 1
         Do ip2 = 0, plotdef%parallelogram%grid(2) - 1
            ip = ip + 1
            t1 = dble (ip1) / dble (plotdef%parallelogram%grid(1))
            t2 = dble (ip2) / dble (plotdef%parallelogram%grid(2))
            vpl (:, ip) = t1 * vl1 (:) + t2 * vl2 (:) + vclp2d (:, 1)
         End Do
      End Do
! evaluate the functions at the grid points
      Do i = 1, nf
         Call rfarray (lmax, ld, rfmt(:, :, :, i), rfir(:, i), ip, vpl, &
        & fp(:, i))
      End Do
! write the functions to file



      Call xml_NewElement (xf, "title")
      Call xml_AddCharacters (xf, trim(input%title))
      Call xml_endElement (xf, "title")
      Call xml_NewElement (xf, "grid")
      Write (buffer, '(2I6)') plotdef%parallelogram%grid(:)
      Call xml_AddAttribute (xf, "gridticks", trim(adjustl(buffer)))
      write(buffer, '(6G18.10)') plotdef%parallelogram%origin%coord
      call xml_AddAttribute (xf, "origin",trim(adjustl(buffer)))
      !write x axis description
      Call xml_NewElement (xf, "axis")
      call xml_AddAttribute (xf, "name", "a")
        Call xml_AddAttribute (xf, "label", get_label(labels,1))
     Call xml_AddAttribute (xf, "latexunit", get_latexunit(labels,1))
     Call xml_AddAttribute (xf, "graceunit", get_graceunit(labels,1))
      write(buffer, '(6G18.10)') plotdef%parallelogram%pointarray(1)%point%coord
      call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
      delta=(plotdef%parallelogram%pointarray(1)%point%coord&
      &-plotdef%parallelogram%origin%coord)&
      &/ plotdef%parallelogram%grid(1)
      write(buffer, '(6G18.10)')  delta
      call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))
      write(buffer, '(6G18.10)')  DNRM2(3,delta,1)
      call xml_AddAttribute (xf, "deltas", trim(adjustl(buffer)))


      Call xml_endElement (xf, "axis")
      !write y axis description
      Call xml_NewElement (xf, "axis")
      call xml_AddAttribute (xf, "name", "b")
      Call xml_AddAttribute (xf, "label", get_label(labels,2))
      Call xml_AddAttribute (xf, "latexunit", get_latexunit(labels,2))
  Call xml_AddAttribute (xf, "graceunit", get_graceunit(labels,2))
      write(buffer, '(6G18.10)') plotdef%parallelogram%pointarray(2)%point%coord
      call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
      delta=(plotdef%parallelogram%pointarray(2)%point%coord&
      &-plotdef%parallelogram%origin%coord)&
      &/ plotdef%parallelogram%grid(2)
       call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))
       write(buffer, '(6G18.10)')  DNRM2(3,delta,1)
       call xml_AddAttribute (xf, "deltas", trim(adjustl(buffer)))
       Call xml_endElement (xf, "axis")
      call xml_NewElement(xf,"value")
       Call xml_AddAttribute (xf, "label", get_label(labels,3))
       Call xml_AddAttribute (xf, "latexunit", get_latexunit(labels,3))
  Call xml_AddAttribute (xf, "graceunit", get_graceunit(labels,3))
           Call xml_endElement (xf, "value")
       Call xml_endElement (xf, "grid")
       ip=0
       Do ifunction = 1, nf
         Call xml_NewElement (xf, "function")
         call xml_AddAttribute (xf, "name", "")
         Do ip1 = 0, plotdef%parallelogram%grid(1) - 1
            Call xml_NewElement (xf, "row")
            call xml_AddAttribute (xf, "const", "x")
            write(buffer20, '(I14)')  ip1
            call xml_AddAttribute (xf, "index", trim(adjustl(buffer20)))
            Do ip2 = 0, plotdef%parallelogram%grid(2) - 1
               ip=ip+1
               write(buffer20, '(6G18.10)')  fp (ip, ifunction)
               call xml_AddCharacters(xf,buffer20)
            end do
            Call xml_endElement (xf, "row")
         end do
         Call xml_endElement (xf, "function")
      end do

      
      Call xml_Close (xf)
      Deallocate (vpl, fp)
      endif
      Return
End Subroutine
!EOC

