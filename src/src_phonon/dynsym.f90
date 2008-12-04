
! Copyright (C) 2006-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynsym(vpl,dynp)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(inout) :: dynp(3*natmtot,3*natmtot)
! local variables
integer iv(3),isym,lspl,i,j,n
real(8) v1(3),v2(3),s(3,3),t1
! automatic arrays
complex(8) dyns(3*natmtot,3*natmtot)
! external functions
real(8) r3taxi
external r3taxi
! map input vector to first Brillouin zone
v1(:)=vpl(:)
call vecfbz(epslat,bvec,v1,iv)
n=0
dyns(:,:)=0.d0
! use the symmetries which leave vpl invariant
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,v1,v2)
  call vecfbz(epslat,bvec,v2,iv)
  if (r3taxi(v1,v2).lt.epslat) then
    call dynsymapp(isym,v1,dynp,dyns)
    n=n+1
  end if
end do
if (n.eq.0) then
  write(*,*)
  write(*,'("Error(dynsym): no symmetries leave vpl invariant")')
  write(*,*)
  stop
end if
t1=1.d0/dble(n)
dynp(:,:)=t1*dyns(:,:)
! make the matrix Hermitian
do i=1,3*natmtot
  do j=i,3*natmtot
    dynp(i,j)=0.5d0*(dynp(i,j)+conjg(dynp(j,i)))
    dynp(j,i)=conjg(dynp(i,j))
  end do
end do
return
end subroutine

