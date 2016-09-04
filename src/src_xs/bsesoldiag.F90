! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bsesoldiag
! !INTERFACE:
subroutine bsesoldiag(hamsiz, ham, eval, evec)
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   integer :: hamsiz ! Dimension of the hermitian matrix
!   complex(hamsiz,hamsiz) :: ! Upper triangular part of an hermitian matrix
! OUT:
!   real(8) :: eval(hamsiz) ! Real valued eigenvalues in ascending order
!   complex(8) :: evec(hamsiz, hamsiz) ! Corresponding eigenvectors
! !DESCRIPTION:
!   Takes a upper triangular part of an hermitian matrix and finds
!   all eigenvalues and eigenvectors using the lapack routine {\tt zheevr}.

  implicit none

  ! Arguments
  integer, intent(in) :: hamsiz
  complex(8), intent(in) :: ham(hamsiz, hamsiz)
  real(8), intent(out) :: eval(hamsiz)
  complex(8), intent(out) :: evec(hamsiz, hamsiz)

  ! Local variables
  real(8) :: vl, vu, abstol
  integer :: il, iu, neval, lwork, info, lrwork, liwork

  ! Allocatable arrays
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:), ifail(:),isuppz(:)

  ! External functions
  real(8), external :: dlamch

  ! Smallest and largest eigenvalue indices
  ! (not referenced since range is set to all in zheevr)
  il = 1
  iu = hamsiz

  ! Tolerance parameter
  abstol = 2.d0 * dlamch('s')

  ! Workspace size (*** improve later ***)
  lwork = (32+1) * hamsiz

  lrwork=-1 !hamsiz
  liwork=-1 !5*hamsiz
  lwork=-1
  iu=hamsiz
  allocate(work(1), rwork(1), iwork(1), isuppz(1))

  call zheevr('v', 'a', 'u', hamsiz, ham, hamsiz, vl, vu, il, iu, &
    & abstol, neval, eval, evec, hamsiz, isuppz, work, lwork, rwork,&
    & lrwork, iwork, liwork, info)

  ! Adjust workspace
  lrwork=int(rwork(1))
  liwork=int(iwork(1))
  lwork=int(work(1))
  deallocate(work, rwork, iwork, isuppz)
  allocate(work(lwork), rwork(lrwork), iwork(liwork))
  allocate(isuppz(hamsiz*2))

  call zheevr('v', 'a', 'u', hamsiz, ham, hamsiz, vl, vu, il, iu, &
    & abstol, neval, eval, evec, hamsiz, isuppz, work, lwork, rwork,&
    & lrwork, iwork, liwork, info)

  deallocate(isuppz)

  if(info .ne. 0) then
   write(*,*)
   write(*, '("Error(bsesoldiag): zheevx returned non-zero info:", i6)') info
   write(*,*)
   call terminate
  end if

  deallocate(work, rwork, iwork)

end subroutine bsesoldiag

