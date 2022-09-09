subroutine comparefiles(file1,file2,diff,identical)
implicit none
character(256),intent(in)::file1,file2,diff
character s
logical, intent(out):: identical
integer::iostat
character(256)::command
logical existfile1,existfile2,exist
inquire(file=trim(file1),exist=exist)
existfile1=exist
inquire(file=trim(file2),exist=exist)
existfile2=exist
if (existfile1.and.existfile2) then

write(command,*)"diff -b ", trim(file1)," ",trim(file2)," >",trim(diff)
CALL SYSTEM(COMMAND)
open(782,file=trim(diff),action='READ',status='OLD',form='FORMATTED', &
 iostat=iostat)
 if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(comparefiles): error opening ",A)') trim(diff)
  write(*,*)
  stop
end if
read(782,*,END=10)s
	go to 11
	write(*,*)"char", s ,"?"
10 	continue
	close (782,STATUS="DELETE")
	identical=.true.
	go to 12
11 	continue
	identical=.false.
	close(782)
12 	continue
else
identical=.false.
if (.not.existfile1) write(*,*)"There is no ", trim(file1), "to compare"
if (.not.existfile2) write(*,*)"There is no " ,trim(file2), "to compare"
endif
	return
end subroutine
