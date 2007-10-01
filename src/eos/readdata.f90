subroutine readdata
use modeos
implicit none
! local variables
integer ipt
open(50,file='eos.in',action='READ',status='OLD',form='FORMATTED')
read(50,*) cname
read(50,*) natoms
read(50,*) etype
read(50,*) vplt1,vplt2,nvplt
read(50,*) nevpt
allocate(vpt(nevpt),ept(nevpt))
do ipt=1,nevpt
  read(50,*) vpt(ipt),ept(ipt)
end do
close(50)
return
end subroutine
