
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writegeomxml
use modspacegroup
use modsymmetries
use FoX_wxml
implicit none
! local variables
integer::is, ia, ip, i
character(128)::buffer
 type(xmlf_t), save::xf
call xml_OpenFile ("geometry.out.xml", xf, replace=.true., pretty_print=.true.)
call xml_NewElement(xf, "input")
call xml_NewElement(xf, "title")
call xml_addCharacters(xf, "title")
call xml_endElement(xf, "title")
call xml_NewElement(xf, "structure")
if(symmetries%lattice%primcell) buffer="true"
if (.not. symmetries%lattice%primcell)buffer="false"
call xml_AddAttribute(xf, "primcell", trim(adjustl(buffer)))
call xml_AddAttribute(xf, "speciespath", trim(adjustl(symmetries%lattice%speciespath)))
call xml_NewElement(xf, "symmetries")
call xml_AddAttribute(xf, "HermannMauguinSymbol", trim( symmetries%HermannMauguinSymbol))
call xml_AddAttribute(xf, "HallSymbol", trim(hall))
call xml_AddAttribute(xf, "SchoenfliesSymbol", trim(schn))
call xml_AddAttribute(xf, "spaceGroupNumber", trim(num))
call xml_NewElement(xf, "lattice")
write(buffer, '(G18.10)') symmetries%lattice%a
call xml_AddAttribute(xf, "a", trim(adjustl(buffer)))
write(buffer, '(G18.10)') symmetries%lattice%b
call xml_AddAttribute(xf, "b", trim(adjustl(buffer)))
write(buffer, '(G18.10)') symmetries%lattice%c
call xml_AddAttribute(xf, "c", trim(adjustl(buffer)))
write(buffer, '(G18.10)') symmetries%lattice%ab
call xml_AddAttribute(xf, "ab", trim(adjustl(buffer)))
write(buffer, '(G18.10)') symmetries%lattice%ac
call xml_AddAttribute(xf, "ac", trim(adjustl(buffer)))
write(buffer, '(G18.10)') symmetries%lattice%bc
call xml_AddAttribute(xf, "bc", trim(adjustl(buffer)))
write(buffer, *) symmetries%lattice%ncell
call xml_AddAttribute(xf, "ncell", trim(adjustl(buffer)))
call xml_endElement(xf, "lattice")
call xml_NewElement(xf, "WyckoffPositions")
do is=1, nspecies
call xml_NewElement(xf, "wspecies")
call xml_AddAttribute(xf, "speciesfile",  trim(adjustl( symmetries%WyckoffPositions%&
wspeciesarray(is)%wspecies%speciesfile)))
  do ip=1, nwpos(is)
     call xml_NewElement(xf, "wpos")
     write(buffer,'(3G18.10)')&
      symmetries%WyckoffPositions%wspeciesarray(is)%wspecies%wposarray(ip)%wpos%coord(:)
     call xml_AddAttribute(xf, "coord", trim(adjustl(buffer)))
     call xml_endElement(xf, "wpos")
  end do
  call xml_endElement(xf, "wspecies")
end do
call xml_endElement(xf, "WyckoffPositions")
call xml_endElement(xf, "symmetries")
call xml_NewElement(xf, "crystal")
if ( symmetries%lattice%scale.ne.1) then
	write(buffer,'(G18.10)') symmetries%lattice%scale
    call xml_AddAttribute(xf, "scale", trim(adjustl(buffer)))
endif
if ( symmetries%lattice%stretch(1).ne.1 .or.&
 symmetries%lattice%stretch(2).ne.1 .or.&
 symmetries%lattice%stretch(3).ne.1 ) then
     write(buffer,'(3G18.10)') symmetries%lattice%stretch
     call xml_AddAttribute(xf, "stretch", trim(adjustl(buffer)))
endif
do i=1, 3
call xml_newElement(xf, "basevect")
write(buffer, '(3G18.10)') avecnew(:, i)
call  xml_AddCharacters(xf, trim(buffer))
call xml_endElement(xf, "basevect")
enddo
call xml_endElement(xf, "crystal")
do is=1, nspecies
call xml_NewElement(xf, "species")
call xml_AddAttribute(xf, "speciesfile", &
 trim(adjustl(   symmetries%WyckoffPositions%wspeciesarray(is)%wspecies%speciesfile)))
  do ia=1, natoms(is)
     call xml_NewElement(xf, "atom")
     write(buffer, '(3G18.10)')atposlnew(:,ia,is)
     call xml_AddAttribute(xf, "coord", trim(adjustl(buffer)))
     call xml_endElement(xf, "atom")
  end do
  call xml_endElement(xf, "species")
end do
call xml_endElement(xf, "structure")
 call xml_NewElement(xf, "groundstate")
  call xml_AddAttribute(xf, "ngridk", " 2 2 2")
write(*, *)
write(*, '("Info(writegeom):")')
write(*, '(" EXCITING lattice vectors and atomic positions written to &
 &geometry.out.xml")')
 call xml_close(xf)
return
end subroutine
