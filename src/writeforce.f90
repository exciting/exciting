
! Copyright (C) 2004-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

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
  write(fnum,'(" species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("  atom : ",I4)') ia
    write(fnum,'("   Hellmann-Feynman",T30,": ",3F14.8)') forcehf(:,ias)
    write(fnum,'("   core correction",T30,": ",3F14.8)') forcecr(:,ias)
    write(fnum,'("   IBS",T30,": ",3F14.8)') forceibs(:,ias)
    write(fnum,'("   total force",T30,": ",3F14.8)') forcetot(:,ias)
    t1=sqrt(forcetot(1,ias)**2+forcetot(2,ias)**2+forcetot(3,ias)**2)
    write(fnum,'("   total magnitude",T30,": ",F14.8)') t1
  end do
end do
call flushifc(fnum)
return
end subroutine

