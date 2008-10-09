program test
use modreport
use modmpi
testplan_name ="test2"
!call inittestoutputfile()

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
!call system("rm *.OUT")
call readinput()
call gndstate()
if (rank.eq.0) then
call write_evec_formatted()
endif
!call finalizeoutput()
end program
