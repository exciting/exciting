subroutine writeforce(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
real(8) t1
write(fnum,*)
write(fnum,'("Forces :")')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'(" species ",I4," atom ",I4," : ")') is,ia
    write(fnum,'("  Hellmann-Feynman : ",3F14.8)') forcehf(:,ias)
    write(fnum,'("  core correction  : ",3F14.8)') forcecr(:,ias)
    write(fnum,'("  IBS              : ",3F14.8)') forceibs(:,ias)
    write(fnum,'("  total force      : ",3F14.8)') forcetot(:,ias)
    t1=sqrt(forcetot(1,ias)**2+forcetot(2,ias)**2+forcetot(3,ias)**2)
    write(fnum,'("  total magnitude  : ",F14.8)') t1
  end do
end do
call flushifc(fnum)
return
end subroutine

