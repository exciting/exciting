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

write(reportfileu,*)testunitname,inputf,outputf,status
write(*,*)testunitname,inputf,outputf,status

end subroutine
