subroutine inittestoutputfile
#include "../../../src/version.inc"
use modmain,only: version
use modreport
implicit none
character(256) tmp
write (tmp,*)"../",trim(adjustl (testplan_name)),".xml"!,".",GITHASH
reportfilename="../report.xml"
reportfileu=84
write(*,*)"Result is written to: ",reportfilename
open(reportfileu,file=reportfilename,action="WRITE",form='FORMATTED')
write(reportfileu,'("<?xml version=""1.0"" encoding=""UTF-8""?>")')
write(reportfileu,*)"<report git_commit="" ", GITHASH,""" >"
end subroutine inittestoutputfile
