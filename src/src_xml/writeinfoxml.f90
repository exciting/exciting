subroutine writeinfoxml
use scl_xml_out_Module
use modmain

   call setAttribute(root, "date", trim(adjustl(buffer)))
   call setAttribute(root, "time", trim(adjustl(buffer)))
   call setAttribute(  ngroundstate, "Volume", trim(adjustl(buffer)))
   ncrystal => createElementNS(sclDoc, "", "crystal")
   dummy => appendChild(ngroundstate, crystal)
   do I=1,3
    nbasevec => createElementNS(sclDoc, "", "basevect")
    content=>createTextNode(scldoc, trim(adjustl(buffer)) )
    dummy => appendChild(basevect, content)
   dummy => appendChild(crystal, cbasevect)
   end do
     do I=1,3
    nreziprvec => createElementNS(sclDoc, "", "reziprvect")
    content=>createTextNode(scldoc, trim(adjustl(buffer)) )
    dummy => appendChild(nreciprvec, content)
   dummy => appendChild(ncrystal, nreciprvec)
   end do
   call setAttribute(ncrystal, "unitCellVolume", trim(adjustl(buffer)))
   call setAttribute(ncrystal, "BrillouinZoneVolume", trim(adjustl(buffer)))
end subroutine

