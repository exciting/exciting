program test

call inittestoutputfile()

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
call system("rm *.OUT")

call test_readinput()
call test_gndstate_init()
call test_comparefiles()
call test_GEOMETRY()
call test_EFERMI()
call test_EIGVAL()
call test_EVALCORE()
call test_EQATOMS()
call test_LINENGY()

call finalizeoutput()

end program
