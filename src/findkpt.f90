
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findkpt(vpl,isym,ik)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(out) :: isym
integer, intent(out) :: ik
! local variables
integer lspl,iv(3)
real(8) s(3,3),v1(3),v2(3),t1
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,vpl,v1)
  call r3frac(epslat,v1,iv)
  do ik=1,nkpt
    v2(:)=vkl(:,ik)
    call r3frac(epslat,v2,iv)
    t1=abs(v1(1)-v2(1))+abs(v1(2)-v2(2))+abs(v1(3)-v2(3))
    if (t1.lt.epslat) return
  end do
end do
write(*,*)
write(*,'("Error(findkpt): equivalent k-point not in set")')
write(*,'(" Requested k-point : ",3G18.10)') vpl
write(*,*)
stop
end subroutine

