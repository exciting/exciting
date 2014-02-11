!
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: plot3d
! !INTERFACE:
!
Subroutine plot3d (plotlabels3d, nf, lmax, ld, rfmt, rfir, plotdef)
! !USES:
    use modplotlabels
    use modinput
    use mod_muffin_tin
    use mod_atoms
    use mod_Gvector
    use FoX_wxml
    use modmpi
!
! !INPUT/OUTPUT PARAMETERS:
!   plotlabels : plot file number (in,integer)
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
!   Modified, February 2011 (DIN)
!   Bug fixing, May 2012 (DIN)
!EOP
!BOC
    Implicit None
! arguments
    type(plotlabels), Intent (In) :: plotlabels3d
    Integer, Intent (In) :: nf
    Integer, Intent (In) :: lmax
    Integer, Intent (In) :: ld
    Real (8), Intent (In) :: rfmt (ld, nrmtmax, natmtot, nf)
    Real (8), Intent (In) :: rfir (ngrtot, nf)
    Type (plot3d_type), Intent (In) :: plotdef
! local variables
    Integer :: np, ip1, ip2, ip3, i, ifunction 
    Real (8) :: tmpv(3)
    Character (512) :: buffer
    Character (20) :: buffer20
    Type (xmlf_t), Save :: xf
! allocatable arrays
    Integer,  Allocatable :: ipmap (:,:,:)
    Real (8), Allocatable :: vpl (:,:)
    Real (8), Allocatable :: fp (:, :)

!-------------------------------------------------------------------------------

    If (rank .Eq. 0) Then
      
        If ((nf .Lt. 1) .Or. (nf .Gt. 4)) Then
            Write(*,*)
            Write(*,'("Error(plot3d): invalid number of functions : ", I8)') nf
            Write(*,*)
            Stop
        End If
!
! allocate the grid point arrays
!
        allocate(ipmap(0:plotdef%box%grid(1),&
        &  0:plotdef%box%grid(2), 0:plotdef%box%grid(3)))
        allocate(vpl(3,(plotdef%box%grid(1)+1)* &
        &  (plotdef%box%grid(2)+1)*(plotdef%box%grid(3)+1)))
!
! generate the 3d point grid and reduce it using the crystal symmetry
!
        call gengrid(plotdef%box%grid,np,ipmap,vpl)
!      
! evaluate the total density at the reduced grid points
!
        allocate (fp(np,nf))
        do i = 1, nf
            call rfarray(lmax, ld, rfmt(:, :, :, i), rfir(:, i), np, vpl, fp(:, i))
        end do
!
! write xml
!
        write (buffer,*) plotlabels3d%filename,".xml"
        call xml_OpenFile( adjustl(trim(buffer)) , xf, replace=.True., pretty_print=.True.)
        call xml_NewElement(xf, "plot3d")
        call xml_NewElement(xf, "title")
        call xml_AddCharacters(xf, trim(input%title))
        call xml_endElement(xf, "title")
        call xml_NewElement(xf, "grid")
        write (buffer,'(3I6)') plotdef%box%grid(1)+1, &
       &  plotdef%box%grid(2)+1, plotdef%box%grid(3)+1
        call xml_AddAttribute(xf, "gridticks", trim(adjustl(buffer)))
        write (buffer,'(3F12.3)') plotdef%box%origin%coord
        call xml_AddAttribute(xf, "origin", trim(adjustl(buffer)))
! write x axis description
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  plotdef%box%origin%coord,1,0.d0,tmpv(1),1)
        write(buffer,'(3F12.3)') tmpv
        call xml_AddAttribute(xf, "originrs", trim(adjustl(buffer)))
        call xml_NewElement(xf, "axis")
        call xml_AddAttribute(xf, "name", "a")
        call xml_AddAttribute(xf, "label", get_label(plotlabels3d,1))
        call xml_AddAttribute (xf, "latexunit", get_latexunit(plotlabels3d,1))
        call xml_AddAttribute (xf, "graceunit", get_graceunit(plotlabels3d,1))
        write(buffer,'(3F12.3)') plotdef%box%pointarray(1)%point%coord
        call xml_AddAttribute(xf, "endpoint", trim(adjustl(buffer)))
        write(buffer,'(3F12.3)') &
       &  (plotdef%box%pointarray(1)%point%coord-plotdef%box%origin%coord)/ &
       &  plotdef%box%grid(1)
        call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  plotdef%box%pointarray(1)%point%coord,1,0.d0,tmpv(1),1)
        write (buffer,'(3F12.3)') tmpv
        call xml_AddAttribute(xf, "endpointrs", trim(adjustl(buffer)))
        call xml_endElement(xf, "axis")
! write y axis description
        call xml_NewElement(xf, "axis")
        call xml_AddAttribute(xf, "name", "b")
        call xml_AddAttribute(xf, "label", get_label(plotlabels3d,2))
        call xml_AddAttribute(xf, "latexunit", get_latexunit(plotlabels3d,2))
        call xml_AddAttribute(xf, "graceunit", get_graceunit(plotlabels3d,2))
        write (buffer,'(3F12.3)') plotdef%box%pointarray(2)%point%coord
        call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
        write (buffer,'(3F12.3)') &
       &  (plotdef%box%pointarray(2)%point%coord-plotdef%box%origin%coord)/ &
       &  plotdef%box%grid(2)
        call xml_AddAttribute(xf, "delta", trim(adjustl(buffer)))
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  plotdef%box%pointarray(2)%point%coord,1,0.d0,tmpv(1),1)
        write(buffer, '(3F12.3)') tmpv
        call xml_AddAttribute(xf, "endpointrs", trim(adjustl(buffer)))
        call xml_endElement(xf, "axis")
! write z axis description
        call xml_NewElement(xf, "axis")
        call xml_AddAttribute(xf, "name", "c")
        call xml_AddAttribute(xf, "label", get_label(plotlabels3d,3))
        call xml_AddAttribute(xf, "latexunit", get_latexunit(plotlabels3d,3))
        call xml_AddAttribute(xf, "graceunit", get_graceunit(plotlabels3d,3))
        write(buffer,'(3F12.3)') plotdef%box%pointarray(3)%point%coord
        call xml_AddAttribute(xf, "endpoint", trim(adjustl(buffer)))
        write(buffer, '(3F12.3)') &
       &  (plotdef%box%pointarray(3)%point%coord-plotdef%box%origin%coord)/ &
       &  plotdef%box%grid(3)
        call xml_AddAttribute(xf, "delta", trim(adjustl(buffer)))
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  plotdef%box%pointarray(3)%point%coord,1,0.d0,tmpv(1),1)
        write (buffer, '(3F12.3)') tmpv
        call xml_AddAttribute (xf, "endpointrs", trim(adjustl(buffer)))
        call xml_endElement(xf, "axis")
!
! write functions to file
!        
        call xml_NewElement(xf,"value")
        call xml_AddAttribute(xf, "label", get_label(plotlabels3d,4))
        call xml_AddAttribute(xf, "latexunit", get_latexunit(plotlabels3d,4))
        call xml_AddAttribute(xf, "graceunit", get_graceunit(plotlabels3d,4))
        call xml_endElement(xf,"value")
        call xml_endElement(xf, "grid")
        do ifunction = 1, nf
            call xml_NewElement(xf, "function")
            write(buffer20,'(I14)') np
            call xml_AddAttribute(xf, "n", trim(adjustl(buffer20)))
            do ip3 = 0, plotdef%box%grid(3)
                call xml_NewElement(xf, "row")
                call xml_AddAttribute(xf, "const", "c")
                write(buffer20,'(I14)') ip3
                call xml_AddAttribute(xf, "index", trim(adjustl(buffer20)))
                do ip2 = 0, plotdef%box%grid(2)
                    call xml_NewElement(xf, "row")
                    call xml_AddAttribute(xf, "const", "b")
                    write(buffer20,'(I14)') ip2
                    call xml_AddAttribute(xf, "index", trim(adjustl(buffer20)))
                    do ip1 = 0, plotdef%box%grid(1)
                        write(buffer20,'(6G18.10)') fp(ipmap(ip1,ip2,ip3),ifunction)
                        call xml_AddCharacters(xf, buffer20)
                    end do
                    call xml_endElement(xf, "row")
                end do
                call xml_endElement(xf, "row")
            end do
            call xml_NewElement(xf, "function")
        end do

        deallocate(vpl, fp)
        deallocate(ipmap)
        call xml_Close(xf)
      
    endif
    
CONTAINS
!
!==================================================================
!
Subroutine gengrid (ngridp, npt, ipmap, vpl)
!
! Created, February 2011 (D. Nabok) 
! 3D real space grid is reduced using the system symmetry
!
      Use modinput
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In)  :: ngridp(3)
      Integer, Intent (Out) :: npt
      Integer, Intent (Out) :: ipmap(0:ngridp(1), 0:ngridp(2), 0:ngridp(3))
      Real (8), Intent (Out) :: vpl(3,(ngridp(1)+1)*(ngridp(2)+1)*(ngridp(3)+1))
! local variables
      Integer :: i1, i2, i3, ip, jp
      Integer :: isym, lspl
      integer :: iv(3)
      Real (8) :: v1(3), v2(3), v3(3)
      Real (8) :: s(3,3), t1

      !-----------------------------------------------------------------

      ip = 0
      Do i3 = 0, ngridp(3)
         v1(3) = dble(i3)/dble(ngridp(3))
         Do i2 = 0, ngridp(2)
            v1(2) = dble(i2)/dble(ngridp(2))
            Do i1 = 0, ngridp(1)
               v1(1) = dble(i1)/dble(ngridp(1))
               ! determine if this point is equivalent to one already in the set
               Do isym = 1, nsymcrys
                  v2(:) = v1(:)+vtlsymc(:,isym)+plotdef%box%origin%coord(:)
                  Call r3frac(input%structure%epslat,v2,iv)
                  lspl = lsplsymc (isym)
                  s(:,:) = dble(symlatc(:,:,lspl))
                  Call r3mv(s,v2,v3)
                  Call r3frac(input%structure%epslat,v3,iv)
                  Do jp = 1, ip
                     t1 = Abs(vpl(1,jp)-v3(1)) + &
                    &     Abs(vpl(2,jp)-v3(2)) + &
                    &     Abs(vpl(3,jp)-v3(3))
                     If (t1 .Lt. input%structure%epslat) Then
                        ! equivalent point found
                        ipmap(i1,i2,i3) = jp
                        goto 10
                     End If
                  End Do
               End Do
               ! add new point to set
               ip = ip+1
               ipmap(i1,i2,i3) = ip
               vpl(:,ip) = v2(:)
10             Continue
            End Do
         End Do
      End Do
      npt = ip
      Return
End Subroutine gengrid


End Subroutine
!EOC
!
