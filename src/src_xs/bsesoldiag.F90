
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bsesoldiag(hamsiz,nev,ham,eval,evec)
  implicit none
  ! arguments
  integer, intent(in) :: hamsiz,nev
  complex(8), intent(in) :: ham(hamsiz,hamsiz)
  real(8), intent(out) :: eval(nev)
  complex(8), intent(out) :: evec(hamsiz,nev)
  ! local variables
  real(8) :: vl,vu,abstol
  integer :: il,iu,neval,lwork,info
  ! allocatable arrays
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:),ifail(:)
  ! external functions
  real(8), external :: dlamch
  ! smallest and largest eigenvalue indices
  il=1
  iu=nev
  ! tolerance parameter
  abstol=2.d0*dlamch('S')
  ! workspace size (*** improve later ***)
  lwork=(32+1)*hamsiz
  allocate(work(lwork),rwork(7*hamsiz),iwork(5*hamsiz),ifail(hamsiz))
  ! LAPACK 3.0 call
  call zheevx('V','I','U',hamsiz,ham,hamsiz,vl,vu,il,iu,abstol,neval,eval, &
       evec,hamsiz,work,lwork,rwork,iwork,ifail,info)
  if (info.ne.0) then
     write(*,*)
     write(*,'("Error(bsesoldiag): zheevx returned non-zero info:",i6)') info
     write(*,*)
     call terminate
  end if
  ! deallocate Hamiltonian array
  deallocate(work,rwork,iwork,ifail)
end subroutine bsesoldiag
