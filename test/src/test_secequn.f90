subroutine test_secequn()
use modreport
implicit none
testunitname="secequn"
inputf="system.in"
outputf="eigenval.out"
call testreport(.false.)

end subroutine

