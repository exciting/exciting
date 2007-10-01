subroutine readinput
use modmain
implicit none
! local variables
integer is,ip
open(50,file='spacegroup.in',action='READ',status='OLD',form='FORMATTED')
! read the Hermann-Mauguin symbol
read(50,*) hrmg
hrmg=adjustl(hrmg)
! read lattice parameters
read(50,*) a,b,c
read(50,*) ab,ac,bc
! read number of unit cells
read(50,*) ncell
if ((ncell(1).lt.1).or.(ncell(2).lt.1).or.(ncell(3).lt.1)) then
  write(*,*)
  write(*,'("Error(readinput): invalid ncell : ",3I8)') ncell
  write(*,*)
  stop
end if
read(50,*) primcell
read(50,*) nspecies
if (nspecies.le.0) then
  write(*,*)
  write(*,'("Error(readinput): nspecies <= 0 : ",I8)') nspecies
  write(*,*)
  stop
end if
if (nspecies.gt.maxspecies) then
  write(*,*)
  write(*,'("Error(readinput): nspecies too large : ",I8)') nspecies
  write(*,'("Adjust maxspecies and recompile code")')
  write(*,*)
  stop
end if
do is=1,nspecies
  read(50,*) spsymb(is),spfname(is)
  read(50,*) nwpos(is)
  if (nwpos(is).le.0) then
    write(*,*)
    write(*,'("Error(readinput): nwpos <= 0 : ",I8)') nwpos(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (nwpos(is).gt.maxwpos) then
    write(*,*)
    write(*,'("Error(readinput): nwpos too large : ",I8)') nwpos(is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxwpos and reompile code")')
    write(*,*)
    stop
  end if
  do ip=1,nwpos(is)
    read(50,*) wpos(:,ip,is)
  end do
end do
close(50)
return
end subroutine

