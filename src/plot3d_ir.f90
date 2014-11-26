!
Subroutine plot3d_ir (plotlabels3d, nf, lmax, ld, rfmt, rfir, plotdef)
! !USES:
    use modplotlabels
    use modinput
    use FoX_wxml
    use modmpi
    use modmain
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
    Integer :: np, ip1, ip2, ip3, i, ifunction, ip, ip_ 
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
    
    v1(:) = (/1.d0, 0.d0, 0.d0/)
    v2(:) = (/0.d0, 1.d0, 0.d0/)
    v3(:) = (/0.d0, 0.d0, 1.d0/)
!
! generate the 3d point grid
!
    np = (ngrid(1)+1)*(ngrid(2)+1)*(ngrid(3)+1)
    allocate(vpl(3,np))
    ip = 0
    Do ip3 = 0, ngrid(3)
       t3 = dble (ip3) / dble (ngrid(3))
       Do ip2 = 0, ngrid(2)
          t2 = dble (ip2) / dble (ngrid(2))
          Do ip1 = 0, ngrid(1)
             t1 = dble (ip1) / dble (ngrid(1))
             ip = ip + 1
             vpl(:, ip) = t1 * v1(:) + &
             &            t2 * v2(:) + &
             &            t3 * v3(:)
          End Do
       End Do
    End Do
!      
! evaluate the total density at the reduced grid points
!

    allocate(fp(np,1))
    ip = 0
    Do ip3 = 0, ngrid(3)
       Do ip2 = 0, ngrid(2)
          Do ip1 = 0, ngrid(1)
             ip  = ip + 1
             ip_ = mod(ip1,ngrid(1))+mod(ip2,ngrid(2))*ngrid(1)+mod(ip3,ngrid(3))*ngrid(1)*ngrid(2)+1
             fp(ip,1) = rfir(ip_,1)*cfunir(ip_)
          End Do
       End Do
    End Do
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
    write (buffer,'(3I6)') ngrid(1)+1, ngrid(2)+1, ngrid(3)+1
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
    write(buffer,'(3F12.3)') v1/ngrid(1)
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
    write(buffer,'(3F12.3)') v2/ngrid(2)
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
    write(buffer,'(3F12.3)') v3/ngrid(3)
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
    call xml_NewElement(xf, "function")
    write(buffer20,'(I14)') np
    call xml_AddAttribute(xf, "n", trim(adjustl(buffer20)))
    ip = 0
    do ip3 = 0, ngrid(3)
        call xml_NewElement(xf, "row")
        call xml_AddAttribute(xf, "const", "c")
        write(buffer20,'(I14)') ip3
        call xml_AddAttribute(xf, "index", trim(adjustl(buffer20)))
        do ip2 = 0, ngrid(2)
            call xml_NewElement(xf, "row")
            call xml_AddAttribute(xf, "const", "b")
            write(buffer20,'(I14)') ip2
            call xml_AddAttribute(xf, "index", trim(adjustl(buffer20)))
            do ip1 = 0, ngrid(1)
                ip = ip+1
                write(buffer20,'(6G18.10)') fp(ip,1)
                call xml_AddCharacters(xf, buffer20)
            end do
            call xml_endElement(xf, "row")
        end do
        call xml_endElement(xf, "row")
    end do
    call xml_endElement(xf, "function")
    deallocate(vpl, fp)
    call xml_Close(xf)

End Subroutine
