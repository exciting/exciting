subroutine updatespecies()

  Use modinput
  Use modmain
  Use modspdeflist
  Use modspdb
  use FoX_wxml
  Use FoX_dom
  
  implicit none
  type(xmlf_t), save :: xf
  character(100) :: buffer
  integer :: is, io, ist, lx, ilo, nlx, i
  
  !---------------------------------------------------------
  ! check first if update of specie files is needed, i.e.,
  ! if the search of linearization energies was used
  !---------------------------------------------------------
  do is = 1, nspecies
    ! apw
    do io = 1, apword(0, is)
      if (apwve(io,0,is)) goto 10
      nlx = size(speziesdeflist(is)%sp%basis%customarray)
      if (nlx>0) then
        do lx = 0, nlx-1
          if (apwve(io,lx,is)) goto 10
        end do
      end if
    end do
    ! lo
    do ilo = 1, nlorb(is)
      do io = 1, lorbord(ilo,is)
        if (lorbve(io,ilo,is)) goto 10
      end do
    end do
  end do
  return
  10 continue
  
  !---------------------------------------------------------
  ! Write down new (updated) xml species file
  !---------------------------------------------------------
  do is = 1, nspecies
  
    call xml_OpenFile (trim(spsymb(is))//'_scf.xml', xf, replace=.true.,pretty_print=.true.)
    call xml_NewElement (xf, "spdb")
    call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XMLSchema-instance", "xsi")
    call xml_AddAttribute(xf, "xsi:noNamespaceSchemaLocation", "../../xml/species.xsd" )
    call xml_NewElement (xf, "sp")
    call xml_AddAttribute(xf, "chemicalSymbol", trim(adjustl(spsymb(is))) )
    call xml_AddAttribute(xf, "name", trim(adjustl(spname(is))) )
    write(buffer,'(G14.6)') spzn(is)
    call xml_AddAttribute(xf, "z", trim(adjustl(buffer)))
    write(buffer,'(G18.10)') spmass(is)
    call xml_AddAttribute(xf, "mass", trim(adjustl(buffer)))
    call xml_NewElement (xf, "muffinTin")
    write(buffer,'(G14.6)') sprmin(is)
    call xml_AddAttribute(xf, "rmin", trim(adjustl(buffer)))
    write(buffer,'(F10.4)') rmt(is)
    call xml_AddAttribute(xf, "radius", trim(adjustl(buffer)))
    write(buffer,'(F10.4)') sprmax(is)
    call xml_AddAttribute(xf, "rinf", trim(adjustl(buffer)))
    write(buffer,*) nrmt(is)
    call xml_AddAttribute(xf, "radialmeshPoints", trim(adjustl(buffer)))
    call xml_endElement(xf, "muffinTin")

    do ist = 1, spnst(is)
      call xml_NewElement (xf, "atomicState")
      write(buffer,*) spn(ist, is) 
      call xml_AddAttribute(xf, "n", trim(adjustl(buffer)))
      write(buffer,*) spl(ist, is)
      call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
      write(buffer,*) spk(ist, is)
      call xml_AddAttribute(xf, "kappa", trim(adjustl(buffer)))
      write(buffer,'(G14.6)') spocc(ist, is)
      call xml_AddAttribute(xf, "occ", trim(adjustl(buffer)))
      if (spcore(ist,is)) then
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

!   Default augmentation
    call xml_NewElement (xf, "default")
    do io = 1, apword(0, is)
      call xml_NewElement (xf, "wf")
      write(buffer,*) apwdm(io, 0, is)
      call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
      write(buffer,'(F8.4)') apwe0(io, 0, is)
      call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
      if (apwve(io, input%groundstate%lmaxapw, is)) then
        call xml_AddAttribute(xf, "searchE", "true")
      else
        call xml_AddAttribute(xf, "searchE", "false")
      endif
      call xml_endElement (xf, "wf")
    end do
    call xml_endElement (xf, "default")

!   Custom augmentation type
    !call xml_AddComment(xf," Custom definitions ")
    nlx = size(speziesdeflist(is)%sp%basis%customarray)
    if (nlx .gt. 0) then
      do lx = 0, nlx-1
        call xml_NewElement (xf, "custom")
        write(buffer,*) lx
        call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
        do io = 1, apword(lx, is)
          call xml_NewElement (xf, "wf")
          write(buffer,*) apwdm(io, lx, is)
          call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
          write(buffer,'(F8.4)') apwe(io, lx, is)
          call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
          if (apwve(io, lx, is)) then
            call xml_AddAttribute(xf, "searchE", "true")
          else
            call xml_AddAttribute(xf, "searchE", "false")
          endif
          call xml_endElement(xf, "wf")
        end do
        call xml_endElement (xf, "custom")
      end do ! lx
    end if

!   Local Orbitals
    !call xml_AddComment(xf," Local Orbitals ")
    do ilo = 1, nlorb(is)
      call xml_NewElement (xf, "lo")
      write(buffer,*) lorbl(ilo, is)
      call xml_AddAttribute(xf, "l", trim(adjustl(buffer)))
      do io = 1, lorbord(ilo, is)
        call xml_NewElement (xf, "wf")
        write(buffer,*) lorbdm(io, ilo, is)
        call xml_AddAttribute(xf, "matchingOrder", trim(adjustl(buffer)))
        write(buffer,'(F8.4)') lorbe(io, ilo, is)
        call xml_AddAttribute(xf, "trialEnergy", trim(adjustl(buffer)))
        if (lorbve(io, ilo, is)) then
          call xml_AddAttribute(xf, "searchE", "true")
        else
          call xml_AddAttribute(xf, "searchE", "false")
        end if
        call xml_endElement (xf, "wf")
      end do ! io
      call xml_endElement (xf, "lo")
    end do ! ilo

    call xml_endElement (xf, "basis")

    call xml_AddComment(xf," This file was automatically generated")
    call xml_AddComment(xf," Default linearization energies have been replaced by that ")
    write(buffer,*) trim(input%groundstate%findlinentype)
    call xml_AddComment(xf," found using findlinentype='"//trim(adjustl(buffer))//"' method ")

    call xml_endElement (xf, "sp")
    call xml_endElement (xf, "spdb")
    call xml_Close(xf)

  end do ! is

end
