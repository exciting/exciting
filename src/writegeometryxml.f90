!BOP
! !ROUTINE: writegeometryxml
! !INTERFACE:

subroutine writegeometryxml(topt)
! !USES:
use modmain
use modinput
use modsp
use  modspdb
use FoX_wxml
! !INPUT/OUTPUT PARAMETERS:
!   topt : if .true. then the filename will be {\ttgeometry_opt.xml}, otherwise
!          {\tt geometry.xml} (in,logical)
! !DESCRIPTION:
!   Outputs the lattice vectors and atomic positions to file, in a format
!   which may be then used directly in {\tt exciting.in}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
logical, intent(in) :: topt
! local variables
integer::is, ia, ip, i
character(128)::buffer
 type(xmlf_t), save::xf
if(topt)then
call xml_OpenFile ("geometry_opt"//trim(filext)//".xml", xf, replace=.true., pretty_print=.true.)
else
call xml_OpenFile ("geometry"//trim(filext)//".xml", xf, replace=.true., pretty_print=.true.)
endif
call xml_NewElement(xf, "input")
call xml_NewElement(xf, "structure")
if(input%structure%primcell) buffer="true"
if (.not. input%structure%primcell)buffer="false"
call xml_AddAttribute(xf, "primcell", trim(adjustl(buffer)))
call xml_AddAttribute(xf, "speciespath", trim(adjustl(input%structure%speciespath)))
call xml_NewElement(xf, "crystal")
do i=1, 3
call xml_newElement(xf, "basevect")
write(buffer, '(3G18.10)') input%structure%crystal%basevect(:,i)
call  xml_AddCharacters(xf, trim(buffer))
call xml_endElement(xf, "basevect")
enddo
call xml_endElement(xf, "crystal")
do is=1, nspecies
call xml_NewElement(xf, "species")
call xml_AddAttribute(xf, "speciesfile", &
&trim(adjustl(input%structure%speciesarray(is)%species%speciesfile)))
 write(buffer,*) int(-1.0 * speziesdeflist(is)%sp%z)
 call xml_AddAttribute(xf, "atomicNumber",trim(adjustl(buffer)))
 call xml_AddAttribute(xf, "chemicalSymbol",trim(adjustl(speziesdeflist(is)%sp%chemicalSymbol)))
  do ia=1, natoms(is)
     call xml_NewElement(xf, "atom")
     write(buffer, '(3G18.10)')input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord
     call xml_AddAttribute(xf, "coord", trim(adjustl(buffer)))
     call xml_endElement(xf, "atom")
  end do
  call xml_endElement(xf, "species")
end do
 call xml_close(xf)
return
end subroutine
