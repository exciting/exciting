
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findkpt(vpl,isym,ik)
use modmain
!<sag>
use modtddft
!</sag>
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(out) :: isym
integer, intent(out) :: ik
! local variables
integer lspl,iv(3)
real(8) s(3,3),v1(3),v2(3)
! external functions
real(8) r3taxi
external r3taxi
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,vpl,v1)
  call r3frac(epslat,v1,iv)
  !<sag>
  if ((task.ge.400).or.(task.le.499)) then
     do ik=1,nkpt0
        v2(:)=vkl0(:,ik)
        call r3frac(epslat,v2,iv)
        if (r3taxi(v1,v2).lt.epslat) return
     end do
  else
  !</sag>
     do ik=1,nkpt
        v2(:)=vkl(:,ik)
        call r3frac(epslat,v2,iv)
        if (r3taxi(v1,v2).lt.epslat) return
     end do
  !<sag>
  end if
  !</sag>
end do
write(*,*)
write(*,'("Error(findkpt): equivalent k-point not in set")')
write(*,'(" Requested k-point : ",3G18.10)') vpl
write(*,*)
stop
end subroutine

