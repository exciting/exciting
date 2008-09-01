program test
use modreport
testplan_name ="test2"
!call inittestoutputfile()

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
!call system("rm *.OUT")
call readinput()
call gndstate()
call write_evec_formatted()

!call finalizeoutput()
end program
