!
!
!
!
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: plot3d
! !INTERFACE:
!
!
Subroutine plot3d (fname, nf, lmax, ld, rfmt, rfir, plotdef)
! !USES:
      Use modinput
      Use modmain
      Use FoX_wxml
      use modmpi

! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngrtot,nf))
! !DESCRIPTION:
!   Produces a 3D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} in the parallelepiped defined by the corner vertices in the
!   global array {\tt vclp3d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Modified, October 2008 (F. Bultmark, F. Cricchio, L. Nordstrom)
!EOP
!BOC
      Implicit None
! arguments
      Character (Len=*), Intent (In) :: fname
      Integer, Intent (In) :: nf
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: rfmt (ld, nrmtmax, natmtot, nf)
      Real (8), Intent (In) :: rfir (ngrtot, nf)
      Type (plot3d_type), Intent (In) :: plotdef
! local variables
      Integer :: np, ip, ip1, ip2, ip3, i, ifunction, fnum = 50
      Real (8) :: v1 (3), v2 (3), v3 (3),tmpv(3)
      Real (8) :: t1, t2, t3
      Character (512) :: buffer, buffer1
      Character (20) :: buffer20
      Type (xmlf_t), Save :: xf
! allocatable arrays
      Real (8), Allocatable :: vpl (:, :)
      Real (8), Allocatable :: fp (:, :)
      buffer = fname // "3d.xml"
 !
!
 If (rank .Eq. 0) Then
      If ((nf .Lt. 1) .Or. (nf .Gt. 4)) Then
         Write (*,*)
         Write (*, '("Error(plot3d): invalid number of functions : ", I&
        &8)') nf
         Write (*,*)
         Stop
      End If
! allocate local arrays
      Allocate (vpl(3, &
     & plotdef%box%grid(1)*plotdef%box%grid(2)*plotdef%box%grid(3)))
      Allocate &
     & (fp(plotdef%box%grid(1)*plotdef%box%grid(2)*plotdef%box%grid(3), &
     & nf))
! generate 3D grid
      v1 (:) = plotdef%box%pointarray(1)%point%coord - &
     & plotdef%box%origin%coord
      v2 (:) = plotdef%box%pointarray(2)%point%coord - &
     & plotdef%box%origin%coord
      v3 (:) = plotdef%box%pointarray(3)%point%coord - &
     & plotdef%box%origin%coord

      ip = 0
      Do ip3 = 0, plotdef%box%grid(3) - 1
         t3 = dble (ip3) / dble (plotdef%box%grid(3))
         Do ip2 = 0, plotdef%box%grid(2) - 1
            t2 = dble (ip2) / dble (plotdef%box%grid(2))
            Do ip1 = 0, plotdef%box%grid(1) - 1
               t1 = dble (ip1) / dble (plotdef%box%grid(1))
               ip = ip + 1
               vpl (:, ip) = t1 * v1 (:) + t2 * v2 (:) + t3 * v3 (:) + &
              & plotdef%box%origin%coord
            End Do
         End Do
      End Do
      np = ip
! evaluate the functions at the grid points
      Do i = 1, nf
         Call rfarray (lmax, ld, rfmt(:, :, :, i), rfir(:, i), np, vpl, &
        & fp(:, i))
      End Do
!write xml
      Call xml_OpenFile (fname//"3d.xml", xf, replace=.True., &
     & pretty_print=.True.)
      Call xml_NewElement (xf, "plot3d")
      Call xml_NewElement (xf, "title")
      Call xml_AddCharacters (xf, trim(input%title))
      Call xml_endElement (xf, "title")
      Call xml_NewElement (xf, "grid")
      Write (buffer, '(3I6)') plotdef%box%grid
      Call xml_AddAttribute (xf, "gridticks", trim(adjustl(buffer)))
      Write (buffer, '(6G18.10)') plotdef%box%origin%coord
      Call xml_AddAttribute (xf, "origin", trim(adjustl(buffer)))
      !write x axis description
       call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3,&
      plotdef%box%origin%coord,1,0.d0,tmpv(1),1)
      Write (buffer, '(3G18.10)') tmpv
      Call xml_AddAttribute (xf, "originrs", trim(adjustl(buffer)))
      Call xml_NewElement (xf, "axis")
      Call xml_AddAttribute (xf, "name", "a")
      Write (buffer, '(6G18.10)') plotdef%box%pointarray(1)%point%coord
      Call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
      Write (buffer, '(6G18.10)') &
     & (plotdef%box%pointarray(1)%point%coord-plotdef%box%origin%coord) &
     & / plotdef%box%grid(1)
      Call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))
      call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3,&
      plotdef%box%pointarray(1)%point%coord,1,0.d0,tmpv(1),1)
      Write (buffer, '(3G18.10)') tmpv
      Call xml_AddAttribute (xf, "endpointrs", trim(adjustl(buffer)))
      Call xml_endElement (xf, "axis")
      !write y axis description
      Call xml_NewElement (xf, "axis")
      Call xml_AddAttribute (xf, "name", "b")
      Write (buffer, '(6G18.10)') plotdef%box%pointarray(2)%point%coord
      Call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
      Write (buffer, '(6G18.10)') &
     & (plotdef%box%pointarray(2)%point%coord-plotdef%box%origin%coord) &
     & / plotdef%box%grid(2)
      Call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))

      call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3,&
      plotdef%box%pointarray(2)%point%coord,1,0.d0,tmpv(1),1)
      Write (buffer, '(3G18.10)') tmpv
      Call xml_AddAttribute (xf, "endpointrs", trim(adjustl(buffer)))
!
      Call xml_endElement (xf, "axis")
            !write z axis description
      Call xml_NewElement (xf, "axis")
      Call xml_AddAttribute (xf, "name", "c")
      Write (buffer, '(6G18.10)') plotdef%box%pointarray(3)%point%coord
      Call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
      Write (buffer, '(6G18.10)') &
     & (plotdef%box%pointarray(3)%point%coord-plotdef%box%origin%coord) &
     & / plotdef%box%grid(3)
      Call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))
       call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3,&
      plotdef%box%pointarray(3)%point%coord,1,0.d0,tmpv(1),1)
      Write (buffer, '(3G18.10)') tmpv
      Call xml_AddAttribute (xf, "endpointrs", trim(adjustl(buffer)))
!
      Call xml_endElement (xf, "axis")
      Call xml_endElement (xf, "grid")
!
!
! write functions to file
      ip = 0
      Do ifunction = 1, nf
         Call xml_NewElement (xf, "function")
         Write (buffer20, '(I14)') np
         Call xml_AddAttribute (xf, "n", trim(adjustl(buffer20)))
         Do ip1 = 0, plotdef%box%grid(3) - 1
            Call xml_NewElement (xf, "row")
            Call xml_AddAttribute (xf, "const", "c")
            Write (buffer20, '(I14)') ip1
            Call xml_AddAttribute (xf, "index", &
           & trim(adjustl(buffer20)))
             Do ip2 = 0, plotdef%box%grid(2) - 1
               Call xml_NewElement (xf, "row")
               Call xml_AddAttribute (xf, "const", "b")
               Write (buffer20, '(I14)') ip2
               Call xml_AddAttribute (xf, "index", &
              & trim(adjustl(buffer20)))
               Do ip3 = 0, plotdef%box%grid(1) - 1
                  ip = ip + 1
                  Write (buffer20, '(6G18.10)') fp (ip, ifunction)
                  Call xml_AddCharacters (xf, buffer20)
               End Do
               Call xml_endElement (xf, "row")
            End Do
            Call xml_endElement (xf, "row")
         End Do
!
         Call xml_NewElement (xf, "function")
      End Do
!
!
      Deallocate (vpl, fp)
      Call xml_Close (xf)
      Close (fnum)
      endif
      Return
End Subroutine
!EOC
