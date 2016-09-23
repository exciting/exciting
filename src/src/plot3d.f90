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
!   Fixed a bug in gengrid, February 2014 (DIN)
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
    Integer :: np, ip1, ip2, ip3, i, ifunction, ip 
    Real (8) :: tmpv(3)
    Character (512) :: buffer
    Character (20) :: buffer20
    Type (xmlf_t), Save :: xf
    Real (8) :: v1(3), v2(3), v3(3), t1, t2, t3
    Real (8) :: boxl(3,4)
! allocatable arrays
    Integer,  Allocatable :: ipmap (:,:,:)
    Real (8), Allocatable :: vpl (:,:)
    Real (8), Allocatable :: fp (:, :)

!-------------------------------------------------------------------------------
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
    v1(:) = plotdef%box%pointarray(1)%point%coord-plotdef%box%origin%coord
    v2(:) = plotdef%box%pointarray(2)%point%coord-plotdef%box%origin%coord
    v3(:) = plotdef%box%pointarray(3)%point%coord-plotdef%box%origin%coord
!
! generate the 3d point grid
!
    If (plotdef%usesym) Then
       boxl(:,1) = plotdef%box%origin%coord
       boxl(:,2) = plotdef%box%pointarray(1)%point%coord-boxl(:,1)
       boxl(:,3) = plotdef%box%pointarray(2)%point%coord-boxl(:,1)
       boxl(:,4) = plotdef%box%pointarray(3)%point%coord-boxl(:,1)
       ! reduce the grid using the crystal symmetry
       call gengrid(plotdef%box%grid,boxl,np,ipmap,vpl)
    Else
       ip = 0
       Do ip3 = 0, plotdef%box%grid(3)
          t3 = dble (ip3) / dble (plotdef%box%grid(3))
          Do ip2 = 0, plotdef%box%grid(2)
             t2 = dble (ip2) / dble (plotdef%box%grid(2))
             Do ip1 = 0, plotdef%box%grid(1)
                t1 = dble (ip1) / dble (plotdef%box%grid(1))
                ip = ip + 1
                ipmap(ip1,ip2,ip3) = ip
                vpl(:, ip) = t1 * v1(:) + &
                &            t2 * v2(:) + &
                &            t3 * v3(:) + &
                &            plotdef%box%origin%coord
             End Do
          End Do
       End Do
       np = ip
    End If
!      
! evaluate the total density at the reduced grid points
!
    allocate(fp(np,nf))
    do i = 1, nf
        call rfarray(lmax, ld, rfmt(:, :, :, i), rfir(:, i), np, vpl, fp(:, i))
    end do
!
! write xml
!
    If (rank .Eq. 0) Then
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
        !write(buffer,'(3F12.3)') plotdef%box%pointarray(1)%point%coord
        write(buffer,'(3F12.3)') v1
        call xml_AddAttribute(xf, "endpoint", trim(adjustl(buffer)))
        write(buffer,'(3F12.3)') v1/plotdef%box%grid(1)
        call xml_AddAttribute (xf, "delta", trim(adjustl(buffer)))
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  v1,1,0.d0,tmpv(1),1)
        write (buffer,'(3F12.3)') tmpv
        call xml_AddAttribute(xf, "endpointrs", trim(adjustl(buffer)))
        call xml_endElement(xf, "axis")
! write y axis description
        call xml_NewElement(xf, "axis")
        call xml_AddAttribute(xf, "name", "b")
        call xml_AddAttribute(xf, "label", get_label(plotlabels3d,2))
        call xml_AddAttribute(xf, "latexunit", get_latexunit(plotlabels3d,2))
        call xml_AddAttribute(xf, "graceunit", get_graceunit(plotlabels3d,2))
        !write (buffer,'(3F12.3)') plotdef%box%pointarray(2)%point%coord
        write(buffer,'(3F12.3)') v2
        call xml_AddAttribute (xf, "endpoint", trim(adjustl(buffer)))
        write(buffer,'(3F12.3)') v2/plotdef%box%grid(2)
        call xml_AddAttribute(xf, "delta", trim(adjustl(buffer)))
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  v2,1,0.d0,tmpv(1),1)
        write(buffer, '(3F12.3)') tmpv
        call xml_AddAttribute(xf, "endpointrs", trim(adjustl(buffer)))
        call xml_endElement(xf, "axis")
! write z axis description
        call xml_NewElement(xf, "axis")
        call xml_AddAttribute(xf, "name", "c")
        call xml_AddAttribute(xf, "label", get_label(plotlabels3d,3))
        call xml_AddAttribute(xf, "latexunit", get_latexunit(plotlabels3d,3))
        call xml_AddAttribute(xf, "graceunit", get_graceunit(plotlabels3d,3))
        !write(buffer,'(3F12.3)') plotdef%box%pointarray(3)%point%coord
        write(buffer,'(3F12.3)') v3
        call xml_AddAttribute(xf, "endpoint", trim(adjustl(buffer)))
        write(buffer,'(3F12.3)') v3/plotdef%box%grid(3)
        call xml_AddAttribute(xf, "delta", trim(adjustl(buffer)))
        call DGEMV('N',3,3,1.d0, input%structure%crystal%basevect(1,1),3, &
       &  v3,1,0.d0,tmpv(1),1)
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
        do i = 1, nf
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
                        write(buffer20,'(6G18.10)') fp(ipmap(ip1,ip2,ip3),i)
                        call xml_AddCharacters(xf, buffer20)
                    end do
                    call xml_endElement(xf, "row")
                end do
                call xml_endElement(xf, "row")
            end do
            call xml_endElement(xf, "function")
        end do

        deallocate(vpl, fp)
        deallocate(ipmap)
        call xml_Close(xf)
      
    end if ! rank==0
    
CONTAINS
!
!==================================================================
!
  Subroutine gengrid (ngridp, b, npt, ipmap, vpl)
!
! Created, February 2011 (D. Nabok) 
! 3D real space grid is reduced using the system symmetry
!
      Use modinput
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In)  :: ngridp(3)
      real(8), intent(in) :: b(3,4)
      Integer, Intent (Out) :: npt
      Integer, Intent (Out) :: ipmap(0:ngridp(1), 0:ngridp(2), 0:ngridp(3))
      Real (8), Intent (Out) :: vpl(3,(ngridp(1)+1)*(ngridp(2)+1)*(ngridp(3)+1))
! local variables
      Integer :: i1, i2, i3, ip, jp
      Integer :: isym, lspl
      integer :: iv(3)
      Real(8) :: r1(3), r2(3), r3(3), vpl_(3)
      Real(8) :: s(3,3), t1

      !-----------------------------------------------------------------
      ip = 0
      Do i3 = 0, ngridp(3)
         r1(3) = dble(i3)/dble(ngridp(3))
         Do i2 = 0, ngridp(2)
            r1(2) = dble(i2)/dble(ngridp(2))
            Do i1 = 0, ngridp(1)
               r1(1) = dble(i1)/dble(ngridp(1))
               ! rescaling
               call r3mv(b(:,2:4),r1,r2)
               ! offset
               r2(:) = r2(:)+b(:,1)
               ! determine if this point is equivalent to one already in the set
               Do isym = 1, nsymcrys
                  lspl = lsplsymc(isym)
                  ! apply symmetry operation S(r+t)
                  s(:,:) = dble(symlat(:,:,lspl))
                  Call r3mv(s,r2(:)+vtlsymc(:,isym),r3)
                  Call r3frac(input%structure%epslat,r3,iv)
                  Do jp = 1, ip
                     vpl_(:) = vpl(:,jp)
                     Call r3frac(input%structure%epslat,vpl_,iv)
                     t1 = Abs(vpl_(1)-r3(1)) + &
                    &     Abs(vpl_(2)-r3(2)) + &
                    &     Abs(vpl_(3)-r3(3))
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
               vpl(:,ip) = r2(:)
10             Continue
            End Do
         End Do
      End Do
      npt = ip
      Return
  End Subroutine gengrid

End Subroutine
!EOC
