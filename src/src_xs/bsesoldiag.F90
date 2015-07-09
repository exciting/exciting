!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine bsesoldiag (hamsiz, nev, ham, eval, evec)
      Implicit None
  ! arguments
      Integer, Intent (In) :: hamsiz, nev
      Complex (8), Intent (In) :: ham (hamsiz, hamsiz)
      Real (8), Intent (Out) :: eval (nev)
      Complex (8), Intent (Out) :: evec (hamsiz, nev)
  ! local variables
      Real (8) :: vl, vu, abstol
      Integer :: il, iu, neval, lwork, info, lrwork,liwork
  ! allocatable arrays
      Complex (8), Allocatable :: work (:)
      Real (8), Allocatable :: rwork (:)
      Integer, Allocatable :: iwork (:), ifail (:),isuppz(:)
  ! external functions
      Real (8), External :: dlamch
  ! smallest and largest eigenvalue indices
      il = 1
      iu = nev
  ! tolerance parameter
      abstol = 2.d0 * dlamch ('S')
  ! workspace size (*** improve later ***)
      lwork = (32+1) * hamsiz
  ! LAPACK 3.0 call
if (.false.) then
      Allocate (work(1), rwork(7*hamsiz), iwork(5*hamsiz), &
     & ifail(hamsiz))
      lwork=-1
      Call zheevx ('V', 'I', 'U', hamsiz, ham, hamsiz, vl, vu, il, iu, &
     & abstol, neval, eval, evec, hamsiz, work, lwork, rwork, iwork, &
     & ifail, info)
      lwork=int(work(1))
      deallocate(work)
      Allocate (work(lwork))

      Call zheevx ('V', 'I', 'U', hamsiz, ham, hamsiz, vl, vu, il, iu, &
     & abstol, neval, eval, evec, hamsiz, work, lwork, rwork, iwork, &
     & ifail, info)
     deallocate(ifail)
else
     lrwork=-1 !hamsiz
     liwork=-1 !5*hamsiz
     lwork=-1
     iu=hamsiz
     Allocate (work(1), rwork(1), iwork(1))
     Call zheevr ('V', 'A', 'U', hamsiz, ham, hamsiz, vl, vu, il, iu, &
     & abstol, neval, eval, evec, hamsiz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, &
     & info)
!     deAllocate (work, rwork, iwork)
     lrwork=int(rwork(1))
     liwork=int(iwork(1))
     lwork=int(work(1))
     deAllocate (work, rwork, iwork)

!     write(*,*) lrwork,liwork,lwork
      Allocate (work(lwork), rwork(lrwork), iwork(liwork))
      allocate(isuppz(hamsiz*2))

      Call zheevr ('V', 'A', 'U', hamsiz, ham, hamsiz, vl, vu, il, iu, &
     & abstol, neval, eval, evec, hamsiz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, &
     & info)
      deallocate(isuppz)
endif

      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(bsesoldiag): zheevx returned non-zero info:&
        &", i6)') info
         Write (*,*)
         Call terminate
      End If
  ! deallocate Hamiltonian array
      Deallocate (work, rwork, iwork)

End Subroutine bsesoldiag

