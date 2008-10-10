subroutine finalizeoutput()
use modreport
implicit none
write(reportfileu,*)"</report>"

close (reportfileu)
end subroutine
