
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findsymeqiv(vpl,vplr,nsc,sc,ivgsc)
  use modmain
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3),vplr(3)
  integer, intent(out) :: nsc,sc(maxsymcrys),ivgsc(3,maxsymcrys)
  ! local variables
  integer :: isym,lspl,iv(3)
  real(8) :: s(3,3),v1(3),t1
  real(8), external :: r3taxi
  ! symmetries that transform non-reduced q-point to reduced one, namely
  ! q1 = s^-1 * q + G_s. Here, q1 is vpl, q is vplr.
  nsc=0
  do isym=1,nsymcrys
     lspl=lsplsymc(isym)
     s(:,:)=dble(symlat(:,:,lspl))
     call r3mtv(s,vplr,v1)
     call r3frac(epslat,v1,iv)
     t1=r3taxi(vpl,v1)
     if (t1.lt.epslat) then
        nsc=nsc+1
        sc(nsc)=isym
        ivgsc(:,nsc)=-iv(:)
     end if
  end do
  if (nsc.eq.0) then
     write(*,*)
     write(*,'("Error(findsymeqiv): p-points are not equivalent by symmetry")')
     write(*,'(" vpl  :",3g18.10)') vpl
     write(*,'(" vplr :",3g18.10)') vplr
     write(*,*)
     call terminate
  end if
end subroutine findsymeqiv
