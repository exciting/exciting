subroutine test_readinput()
use modreport
implicit none
logical passed

testunitname="readinput"
inputf="exciting.in"
outputf="parameters.out"
passed=.false.
call readinput()
passed=.true.
call testreport(passed)

end subroutine

