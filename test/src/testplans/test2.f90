program test
use modreport
testplan_name ="test2"
call inittestoutputfile()

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
call system("rm *.OUT")
call test_readinput()
call test_gndstate_init()
call test_lapacksolver()
call test_arpacksolver()

call finalizeoutput()
end program
