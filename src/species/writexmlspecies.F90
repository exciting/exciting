
! Copyright (C) 2005-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine  writexmlspecies
  use FoX_wxml
  use modspecies
  implicit none
  type(xmlf_t), save::xf
  character(100)::buffer,buffer2
  call xml_OpenFile (trim(spsymb)//'_'//trim(suffix)//'.xml', xf, replace=.true.,pretty_print=.true.)
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
!
! Default augmentation
! 
! basis type
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

! trial energy
  write(buffer,'(F8.4)')boe
  call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))

! search LE
  if(apwve(io)) then
    call xml_AddAttribute(xf, "searchE", "true")
  else
    call xml_AddAttribute(xf, "searchE", "false")
  endif
  call xml_endElement (xf, "default")

!
! Custom augmentation
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
      write(buffer,'(F8.4)') boe
      call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
      if (apwvex(io)) then
        call xml_AddAttribute(xf, "searchE", "true")
      else
        call xml_AddAttribute(xf, "searchE", "false")
      endif
      call xml_endElement (xf, "custom")
    end do ! l
  end if
!
! Local Orbitals
!
   buffer2="false"
   if (fullsearchlocorbsc) buffer2="true"

   if (locorbsc) then
     do ist=1,spnst
       if (.not.spcore(ist)) then
         if ((spl(ist).eq.0).or.(spl(ist).eq.spk(ist))) then
           if (eval(ist).lt.esccut) then
             call xml_NewElement (xf, "lorb")
             write(buffer,*)spl(ist)
             call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
             call xml_NewElement (xf, "wf")
             call xml_AddAttribute(xf, "matchingOrder", "0")
             write(buffer,'(F8.4)')boe
             call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
             call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
             call xml_endElement (xf, "wf")
             call xml_NewElement (xf, "wf")
             call xml_AddAttribute(xf, "matchingOrder", "1")
             write(buffer,'(F8.4)')boe
             call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
             call xml_AddAttribute(xf, "searchE", trim(adjustl(buffer2)))
             call xml_endElement (xf, "wf")
             call xml_NewElement (xf, "wf")
             call xml_AddAttribute(xf, "matchingOrder", "0")
             write(buffer,'(F8.4)') eval(ist)+0.5d0*boe
             call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
             call xml_AddAttribute(xf, "searchE", "true")
             call xml_endElement (xf, "wf")
             call xml_endElement (xf, "lorb")
           end if
         end if
       end if
     end do
   endif ! locorbsc

  call xml_endElement (xf, "basis")

  call xml_AddComment(xf," Exciting code version: "//version//" ")
  call xml_AddComment(xf," Description of method: "//trim(suffix)//" ")
  call xml_AddComment(xf," "//trim(apwdescr)//" ")
  call xml_Close(xf)

end subroutine writexmlspecies
