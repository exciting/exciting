subroutine test_gndstate_init()
use modmain
use modreport
implicit none
logical passed
maxscl=0
passed=.false.
call gndstate()
testunitname="gndstate_init"
inputf="exciting"
outputf="INFO.OUT"
passed=.true.
call testreport(passed)
end subroutine

