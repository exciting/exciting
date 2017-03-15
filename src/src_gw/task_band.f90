
subroutine task_band()

  use modinput
  use modmain
  use modgw
  use modmpi
  Use FoX_wxml

  implicit none

  integer(4)    :: ik, ib, ib0
  real(8)       :: tstart, tend
  integer       :: i, j
  character(80) :: fname, s
  logical       :: exist
  Character (128) :: buffer
  Type (xmlf_t), Save :: xf

  !------------------------
  ! Read KS bandstructure
  !------------------------

  call init0()
  call init1

  fname = 'bandstructure.dat'
  inquire(File=fname, Exist=exist)
  if (.not.exist) then
    write(*,*) 'ERROR(task_band): bandstructure.dat file is not found!'
    write(*,*) '    Run properties/bandstructure first to produce KS spectrum.'
    stop
  end if

  open(70, File='bandstructure.dat', Action='Read', Status='Old')
  read(70,*) s, ib0, nstsv, nkpt
  write(*,*) ib0, nstsv, nkpt
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(dpp1d)) deallocate(dpp1d)
  allocate(dpp1d(nkpt))
  if (allocated(evalsv)) deallocate(evalsv)
  allocate(evalsv(nstsv,nkpt))
  do ib = ib0, nstsv
    do ik = 1, nkpt
      ! Note: evalsv are already shifted to E_f = 0
      read(70,*) i, j, vkl(:,ik), dpp1d(ik), evalsv(ib,ik)
    end do
    read(70,*) ! skip line
  end do
  close(70)

  !--------------------------------------------------------------
  ! read QP energies from file and perform Fourier interpolation 
  !--------------------------------------------------------------
  call getevalqp(nkpt,vkl,evalsv)

  !----------------------------------
  ! write QP bandstructure to disk
  !----------------------------------
  call xml_OpenFile ("bandstructure-qp.xml", xf, replace=.True., &
       & pretty_print=.True.)
  call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
       &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')
  call xml_NewElement( xf, "bandstructure")
  call xml_NewElement( xf, "title")
  call xml_AddCharacters( xf, trim( input%title))
  call xml_endElement( xf, "title")
  open(50, File='BAND-QP.OUT', Action='Write', Form='Formatted')
  open(51, File="bandstructure-qp.dat", Action='Write', Form='Formatted')
  write(51,*) "# ", ibgw, min(nbgw,nstsv), nkpt
  do ib = ibgw, min(nbgw,nstsv)
    call xml_NewElement( xf, "band")
    do ik = 1, nkpt
      ! old format (gwmod-boron) 
      ! write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)
      ! write(51,'(2I6, 5F12.6)') &
      ! &     ib, ik, vkl(:,ik), dpp1d(ik), evalsv(ib,ik)
      ! new format (carbon)
      write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)+efermi
      write(51,'(2I6, 3F12.6, 2G18.10)') &
      &     ib, ik, vkl(:,ik), dpp1d(ik), evalsv(ib,ik)+efermi
      call xml_NewElement( xf, "point")
      write( buffer, '(5G18.10)') dpp1d (ik)
      call xml_AddAttribute( xf, "distance", trim( adjustl( buffer)))
      write( buffer, '(5G18.10)') evalsv( ib, ik)+efermi
      call xml_AddAttribute (xf, "eval", trim( adjustl( buffer)))
      call xml_endElement( xf, "point")
    end do !ik
    call xml_endElement( xf, "band")
    write(50,*)
    write(51,*)
  end do !ib
  close(50)
  close(51)
  call xml_endElement( xf, "bandstructure")
  call xml_close( xf)

  return
end subroutine
