program test
use modreport
use modmain
complex(8), allocatable :: apwalm(:,:,:,:,:)
testplan_name ="test2"
call inittestoutputfile()

! list test routines here and call testreport(testunit,input,output,passed)
! before leaving the routine
call system("rm *.OUT")
call test_readinput()
call test_gndstate_init()
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
          sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))

call test_hamiltonsetup(apwalm)
call test_lapacksolver(apwalm)
call test_arpacksolver(apwalm)
deallocate(apwalm)
call finalizeoutput()
end program
