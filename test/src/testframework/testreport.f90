subroutine testreport(passed)
use modreport
logical,intent(in)::passed

character(6)::status
tests=tests+1
if(passed)then
status="passed"
npassed=npassed+1
else
status="failed"
nfailed=nfailed+1
endif

write(reportfileu,*)"   <test>"
write(reportfileu,*)"        <name>",trim(testunitname),"</name>"
write(reportfileu,*)"        <description>",inputf, outputf,"</description>"
write(reportfileu,*)"        <status>",status,"</status>"
write(reportfileu,*)"        <directory>",trim(tdirectory),"</directory>"
write(reportfileu,*)"   </test>"

write(*,*)testunitname,inputf,outputf,status

end subroutine
