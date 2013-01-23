
! Copyright (C) 2005-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine  writexmlspecies

  use FoX_wxml
  use modspecies

  implicit none
  type(xmlf_t), save::xf
  character(100)::buffer,buffer1,buffer2
  real(8) :: e
  integer :: order

  call xml_OpenFile (trim(spsymb)//trim(suffix)//'.xml', xf, replace=.true.,pretty_print=.true.)
  call xml_NewElement (xf, "spdb")
  call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XMLSchema-instance", "xsi")
  call xml_AddAttribute(xf, "xsi:noNamespaceSchemaLocation", "../../xml/species.xsd" )
  call xml_NewElement (xf, "sp")
  call xml_AddAttribute(xf, "chemicalSymbol", trim(adjustl(spsymb)) )
  call xml_AddAttribute(xf, "name", trim(adjustl(spname)) )
  write(buffer,'(G14.6)')spzn
  call xml_AddAttribute(xf, "z", trim(adjustl(buffer)))
  write(buffer,'(G18.10)')spmass
  call xml_AddAttribute(xf, "mass", trim(adjustl(buffer)))
  call xml_NewElement (xf, "muffinTin")
  write(buffer,'(G14.6)')sprmin
  call xml_AddAttribute(xf, "rmin", trim(adjustl(buffer)))
  write(buffer,'(F10.4)')rmt
  call xml_AddAttribute(xf, "radius", trim(adjustl(buffer)))
  write(buffer,'(F10.4)')sprmax
  call xml_AddAttribute(xf, "rinf", trim(adjustl(buffer)))
  write(buffer,*)nrmt
  call xml_AddAttribute(xf, "radialmeshPoints", trim(adjustl(buffer)))
  call xml_endElement(xf, "muffinTin")

  do ist=1,spnst
     call xml_NewElement (xf, "atomicState")
     write(buffer,*)spn(ist)
     call xml_AddAttribute(xf, "n", trim(adjustl(buffer)))
     write(buffer,*)spl(ist)
     call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
     write(buffer,*)spk(ist)
     call xml_AddAttribute(xf, "kappa", trim(adjustl(buffer)))
     write(buffer,'(G14.6)')spocc(ist)
     call xml_AddAttribute(xf, "occ", trim(adjustl(buffer)))
     if(spcore(ist))then
        buffer="true"
     else
        buffer="false"
     endif
     call xml_AddAttribute(xf, "core", trim(adjustl(buffer)))
     call xml_endElement(xf, "atomicState")
  end do

!--------------------------------------------------------------
! BASIS
!--------------------------------------------------------------
  call xml_NewElement (xf, "basis")

  if (trim(apwdescr).eq.'simple') then
!
!   Default augmentation
! 
!   basis type
    call xml_NewElement (xf, "default")
    if (apword.eq.1) then
      write(buffer,*) 'apw+lo'
    else if (apword.eq.2) then
      write(buffer,*) 'lapw'
    else
      write(*,*) 'Not supported basis type: apword>2'
      stop
    end if
    call xml_AddAttribute(xf, "type", trim(adjustl(buffer)))
!
!   trial energy
    write(buffer,'(F8.4)') boe
    call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
!
!   search LE
    if(apwve(1)) then
      call xml_AddAttribute(xf, "searchE", "true")
    else
      call xml_AddAttribute(xf, "searchE", "false")
    endif
    call xml_endElement (xf, "default")
!
!   Custom augmentation type
!
    if (apwordx.gt.0) then
      do l = 0, maxl
        call xml_NewElement (xf, "custom")
        write(buffer,*) l
        call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
        if (apwordx.eq.1) then
          write(buffer,*) 'apw+lo'
          call xml_AddAttribute(xf, "type", trim(adjustl(buffer)))
        else if (apwordx.eq.2) then
          write(buffer,*) 'lapw'
          call xml_AddAttribute(xf, "type", trim(adjustl(buffer)))
        else
          write(*,*) 'Not supported basis type: apwordx>2'
          stop
        end if
        write(buffer,'(F8.4)') elval(l)
        call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
        if (apwvex(1)) then
          call xml_AddAttribute(xf, "searchE", "true")
        else
          call xml_AddAttribute(xf, "searchE", "false")
        endif
        call xml_endElement (xf, "custom")
      end do ! l
    end if
  
  else if (trim(apwdescr).eq.'advanced') then
!
!   Default augmentation
! 
!   basis type
    call xml_NewElement (xf, "default")
    do io=1,apword
      call xml_NewElement (xf, "wf")
      write(buffer,*) apwdm(io)
      call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
      write(buffer,'(F8.4)') boe
      call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
      if(apwve(io)) then
        call xml_AddAttribute(xf, "searchE", "true")
      else
        call xml_AddAttribute(xf, "searchE", "false")
      endif
      call xml_endElement (xf, "wf")
    end do
    call xml_endElement (xf, "default")
!
!   Custom augmentation type
!
    if (apwordx.gt.0) then
      do l = 0, maxl
        call xml_NewElement (xf, "custom")
        write(buffer,*) l
        call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
        do io = 1, apwordx
          call xml_NewElement (xf, "wf")
          write(buffer,*) apwdmx(io)
          call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
          write(buffer,'(F8.4)') elval(l)
          call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
          if(apwvex(io)) then
            call xml_AddAttribute(xf, "searchE", "true")
          else
            call xml_AddAttribute(xf, "searchE", "false")
          endif
          call xml_endElement(xf, "wf")
        end do
        call xml_endElement (xf, "custom")
!
!       print lo for APW+lo method
!
        if (locorb) then
          buffer2="false"
          if(searchlocorb) buffer2="true"
            call xml_NewElement (xf, "lo")
            write(buffer,*) l
            call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
            call xml_NewElement (xf, "wf")
            call xml_AddAttribute(xf, "matchingOrder", "0")
            write(buffer,'(F8.4)') elval(l)
            call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
            call xml_endElement (xf, "wf")
            call xml_NewElement (xf, "wf")
            call xml_AddAttribute(xf, "matchingOrder", "1")
            write(buffer,'(F8.4)') elval(l)
            call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
            call xml_endElement (xf, "wf")
            call xml_endElement (xf, "lo")
         endif
!
      end do ! l
    end if
  end if
!
! Semi-Core Local Orbitals
!
  if ( locorbsc ) then
    buffer2="false"
    if (searchlocorb) buffer2="true"
    do l = 0, lmax
      do i = 1, nl(l)
        e = el(l,i)
        if (e.lt.esccut) then
            call xml_NewElement (xf, "lo")
            write(buffer,*) l
            call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
            call xml_NewElement (xf, "wf")
            call xml_AddAttribute(xf, "matchingOrder", "0")
            write(buffer,'(F8.4)') elval(l)
            call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
            call xml_endElement (xf, "wf")
            call xml_NewElement (xf, "wf")
            call xml_AddAttribute(xf, "matchingOrder", "1")
            write(buffer,'(F8.4)') elval(l)
            call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
            call xml_endElement (xf, "wf")
            call xml_NewElement (xf, "wf")
            call xml_AddAttribute(xf, "matchingOrder", "0")
            write(buffer,'(F8.4)') e
            call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", "true")
            call xml_endElement (xf, "wf")
            call xml_endElement (xf, "lo")
        end if
      end do
    end do
  endif ! locorbsc
  
!-------------------------------------------------
! Exciting States Local Orbitals
!-------------------------------------------------

  if ( locorbxs ) then
    buffer2="false"
    if (searchlocorbxs) buffer2="true"
    do l = 0, lxsmax
      buffer1="false"
      if ((l<=maxl).and.(apwvex(1))) buffer1="true"
      ! ----
      order=apword
      if (l<=maxl) order=apwordx
      ! ----
      do i = 1, nl(l)
        e = el(l,i)
        if (e.gt.exscut) then
          call xml_NewElement (xf, "lo")
          write(buffer,*) l
          call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
          ! valence radial function
          do io = 0, order-1
            call xml_NewElement (xf, "wf")
            write(buffer,*) io
            call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
            write(buffer,'(F8.4)') elval(l)
            call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
            call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer1)))
            call xml_endElement (xf, "wf")
          end do
          ! xs radial function
          io = 0
          call xml_NewElement (xf, "wf")
          write(buffer,*) io
          call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
          write(buffer,'(F8.4)') e
          call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
          call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
          call xml_endElement (xf, "wf")
          call xml_endElement (xf, "lo")
        end if
      end do ! i
    end do ! l
  endif ! locorbxs

  call xml_endElement (xf, "basis")

  call xml_AddComment(xf," Exciting code version: "//version//" ")
  call xml_AddComment(xf," Generation Scheme: "//trim(apwdescr)//" ")
  call xml_Close(xf)

end subroutine writexmlspecies
