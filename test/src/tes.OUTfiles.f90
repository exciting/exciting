subroutine test_GEOMETRY()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="GEOMETRY.OUT"
file2="../reference/GEOMETRY.OUT"
diff="../output/geomdiff"
call comparefiles(file1,file2,diff,passed)
testunitname="Geometry file"
inputf="GEOMETRY.OUT"
outputf="geomdiff"
call testreport(passed)
end subroutine

subroutine test_EFERMI()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="EFERMI.OUT"
file2="../reference/EFERMI.OUT"
diff="../output/EFERMI.diff"
call comparefiles(file1,file2,diff,passed)
testunitname="EFERMI.OU"
inputf="EFERMI.OUT"
outputf="EFERMI.diff"
call testreport(passed)
end subroutine


subroutine test_EIGVAL()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="EIGVAL.OUT"
file2="../reference/EIGVAL.OUT"
diff="../output/EIGVAL.diff"
call comparefiles(file1,file2,diff,passed)
testunitname="EIGVAL.OU"
inputf="EIGVAL.OUT"
outputf="EIGVAL.diff"
call testreport(passed)
end subroutine


subroutine test_EVALCORE()
use modreport
implicit none
character(256)::file1,file2,diff
logical::passed
file1="EVALCORE.OUT"
file2="../reference/EVALCORE.OUT"
diff="../output/EVALCORE.diff"
call comparefiles(file1,file2,diff,passed)
testunitname="EVALCORE.OU"
inputf="EVALCORE.OUT"
outputf="EVALCORE.diff"
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
testunitname="EQATOMS.OU"
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
testunitname="LINENGY.OU"
inputf="LINENGY.OUT"
outputf="LINENGY.diff"
call testreport(passed)
end subroutine
