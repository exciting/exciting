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
