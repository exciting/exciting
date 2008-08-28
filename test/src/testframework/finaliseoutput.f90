subroutine finalizeoutput()
use modreport
implicit none
write(reportfileu,*)
write(reportfileu,*)
write(reportfileu,*)"-------------------------------------"
write(reportfileu,*) "total tests",tests
write(reportfileu,*) "passed     ",npassed
write(reportfileu,*) "failed     ",nfailed
close (reportfileu)
end subroutine
