subroutine test_GEOMETRY()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="GEOMETRY.OUT"
file2="../reference/GEOMETRY.OUT"
diff="../output/geomdiff"
call comparefiles(file1,file2,diff,passed)
testunitname="Geometry_file"
inputf="GEOMETRY.OUT"
outputf="geomdiff"
call testreport(passed)
end subroutine

subroutine test_EFERMI()
use modreport
implicit none
character(256)::file1,file2,errfile
logical::passed
real:: efermi1,efermi2,tol
tol=1e-8
file1="EFERMI.OUT"
file2="../reference/EFERMI.OUT"
errfile="../output/EFERMI.err"
open(888,file=file1)
read(888,*) efermi1
close(888)
open(888,file=file2)
read(888,*) efermi2
close(888)
if(abs(efermi1-efermi2).lt.tol) then
passed=.true.
else
passed=.false.
endif
open(888,file=errfile)
write(888,*) "error, tol"
write(888,*) efermi1-efermi2, tol
close(888)
testunitname="Fermienergy"
inputf="EFERRMI.OUT"
outputf=errfile
call testreport(passed)
end subroutine




subroutine test_EQATOMS()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="EQATOMS.OUT"
file2="../reference/EQATOMS.OUT"
diff="../output/EQATOMS.diff"
call comparefiles(file1,file2,diff,passed)
testunitname="EQATOMS.OUT"
inputf="EQATOMS.OUT"
outputf="EQATOMS.diff"
call testreport(passed)
end subroutine


subroutine test_LINENGY()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="LINENGY.OUT"
file2="../reference/LINENGY.OUT"
diff="../output/LINENGY.diff"
call comparefiles(file1,file2,diff,passed)
testunitname="LINENGY.OUT"
inputf="LINENGY.OUT"
outputf="LINENGY.diff"
call testreport(passed)
end subroutine
