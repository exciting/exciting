subroutine inittestoutputfile
#include "../../../src/version.inc"
use modmain,only: version
use modreport
implicit none
character(256) tmp
write (tmp,*)"../",trim(adjustl (testplan_name)),".summary"!,".",GITHASH
reportfilename=trim(adjustl(tmp))
reportfileu=84
write(*,*)"Result is written to: ",reportfilename
open(reportfileu,file=reportfilename,action="WRITE",form='FORMATTED')
write(reportfileu,'(" +----------------------------------+")')
write(reportfileu,'(" | EXCITING version ",I1.1,".",I1.1,".",I3.3," started |")') version

write(reportfileu,*)"| git hash id ", GITHASH," |"
write(reportfileu,'(" +----------------------------------+")')

write(reportfileu,*)
write(reportfileu,*)"unit                input               output              status"
write(reportfileu,*)"-------------------------------------------------------------------"
tests=0
nfailed=0
npassed=0
end subroutine inittestoutputfile
