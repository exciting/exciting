program test

call inittestoutputfile(50) !file unit 50

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine

call test_readinput()
call test_gndstate_init()
call test_comparefiles()

call finalizeoutput(50)

end program
