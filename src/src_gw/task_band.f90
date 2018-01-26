
subroutine task_band()

  use modinput
  use modmain
  use modgw
  use modmpi
  Use FoX_wxml

  implicit none

  integer(4)    :: ik, ib, ib0
  real(8)       :: tstart, tend
  integer       :: i, j, is, ia, ias, l
  character(80) :: fname, s
  logical       :: exist, bandchar
  Character (128) :: buffer
  Type (xmlf_t), Save :: xf
  real(8), allocatable :: bc(:,:,:,:)

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
  bandchar = .false.
  if( bandchar) then
    allocate( bc( 0:3, natmtot, nstsv, nkpt))
    bc = 0.d0
    call getevalqp(nkpt,vkl,evalsv,bc)
  else
    call getevalqp(nkpt,vkl,evalsv)
  end if

  !----------------------------------
  ! write QP bandstructure to disk
  !----------------------------------
  call xml_OpenFile ("bandstructure-qp.xml", xf, replace=.True., &
       & pretty_print=.True.)
  call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
       &'/visualizationtemplates/bandstructure2html.xsl" type="text/xsl"')

    if( .not. bandchar) then
      open( 50, file='BAND-QP.OUT', action='WRITE', form='FORMATTED')
      call xml_NewElement( xf, "bandstructure")
      call xml_NewElement( xf, "title")
      call xml_AddCharacters( xf, trim( input%title))
      call xml_endElement( xf, "title")
      open(51, File="bandstructure-qp.dat", Action='Write', Form='Formatted')
      write(51,*) "# ", ibgw, min(nbgw,nstsv), nkpt
      do ib = ibgw, min( nbgw, nstsv)
        call xml_NewElement( xf, "band")
        do ik = 1, nkpt
          ! old format (gwmod-boron) 
          ! write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)
          ! write(51,'(2I6, 5F12.6)') &
          ! &     ib, ik, vkl(:,ik), dpp1d(ik), evalsv(ib,ik)
          ! new format (carbon)
          write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)+efermi
          write(51,'(2I6, 3F12.6, 2G18.10)') ib, ik, vkl(:,ik), dpp1d(ik), evalsv(ib,ik)+efermi
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
      end do
      close(50)
      close(51)
    else
      call xml_NewElement (xf, "bandstructure")
      call xml_AddAttribute (xf, "character", "true")
      call xml_NewElement (xf, "title")
      call xml_AddCharacters (xf, trim(input%title))
      call xml_endElement (xf, "title")
      do is = 1, nspecies
        call xml_NewElement (xf, "species")
        call xml_AddAttribute (xf, "name", trim(spname(is)))
        call xml_AddAttribute (xf, "chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol))
        do ia = 1, natoms (is)
          call xml_NewElement (xf, "atom")
          write (buffer, '(5G18.10)') atposc (:, ia, is)
          call xml_AddAttribute (xf, "coord", &
               & trim(adjustl(buffer)))
          ias = idxas (ia, is)
          write (fname, '("BAND-QP_S", I2.2, "_A", I4.4, ".OUT")') is, ia
          open (50, File=trim(fname), Action='WRITE', Form='FORMATTED')
          !
          do ib = ibgw, min( nbgw, nstsv)
            call xml_NewElement (xf, "band")
            do ik = 1, nkpt
              ! sum band character over l
              call xml_NewElement (xf, "point")
              write (buffer, '(5G18.10)') dpp1d( ik)
              call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
              write (buffer, '(5G18.10)') evalsv( ib, ik)+efermi
              call xml_AddAttribute (xf, "eval", trim(adjustl(buffer)))
              write (buffer, '(5G18.10)') sum( bc( 0:3, ias, ib, ik))
              call xml_AddAttribute (xf, "sum", trim(adjustl(buffer)))
              do l = 0, 3
                call xml_NewElement (xf, "bc")
                write (buffer,*) l
                call xml_AddAttribute (xf, "l", trim(adjustl(buffer)))
                write (buffer, '(5G18.10)') bc( l, ias, ib, ik)
                call xml_AddAttribute (xf, "character", trim(adjustl(buffer)))
                call xml_endElement (xf, "bc")
              end do
              call xml_endElement (xf, "point")
              write (50, '(2G18.10, 20F12.6)') dpp1d( ik), evalsv( ib, ik)+efermi, sum( bc( 0:3, ias, ib, ik)), (bc( l, ias, ib, ik), l=0, 3)
            end do
            call xml_endElement (xf, "band")
            write (50, '("	  ")')
          end do
          call xml_endElement (xf, "atom")
          close (50)
        end do
        call xml_endElement (xf, "species")
      end do
      call xml_endElement( xf, "bandstructure")
      call xml_close( xf)
      deallocate( bc)
    end if  


  return
end subroutine
