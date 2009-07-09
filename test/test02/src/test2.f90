program test
use modreport
use modmpi
use modinput
implicit none
testplan_name ="test2"
!call inittestoutputfile()

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
!call system("rm *.OUT")
 call initMPI

call loadinputDOM()
call setdefault
input=getstructinput(inputnp)
call ifparseerrorstop()
call destroyDOM()
call initatomcounters()
call initlattice
call readspeciesxml
call gndstate()
call finitMPI
if (rank.eq.0) then
call write_evec_formatted()
endif

!call finalizeoutput()
end program
