! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bsesoldiag
! !INTERFACE:
subroutine bsesoldiag(solsize, hamsize, ham, eval, evec)
! !USES:
  use modmpi
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   integer :: solsize ! Number of solutions from lowest EV 
!   integer :: hamsize ! Dimension of the hermitian matrix
!   complex(8) :: ham(hamsize,hamsize) ! Upper triangular part of an hermitian matrix
! OUT:
!   real(8) :: eval(hamsize) ! Real valued eigenvalues in ascending order 
!                            ! (the first solsize elements are set)
!   complex(8) :: evec(hamsize, hamsize) ! Corresponding eigenvectors as columns
!
! !DESCRIPTION:
!   Takes a upper triangular part of an hermitian matrix and finds
!   eigenvalues and eigenvectors using the lapack routine {\tt zheevr}.
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: solsize, hamsize
  complex(8), intent(in) :: ham(hamsize, hamsize)
  real(8), intent(out) :: eval(hamsize)
  complex(8), intent(out) :: evec(hamsize, solsize)

  ! Local variables
  real(8) :: vl, vu, abstol
  integer :: il, iu, lwork, info, lrwork, liwork

  ! Allocatable arrays
  complex(8), allocatable :: work(:)
  real(8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:), isuppz(:)

  ! External functions
  real(8), external :: dlamch

  ! Tolerance parameter
  abstol = 2.d0 * dlamch('s')

  if(solsize < 1 .or. solsize > hamsize) then
    write(*,*) "bsesoldiag (ERROR): Number of requested solutions is &
      & incompatible with size of Hamiltonian."
    write(*,*) "solsize:", solsize
    write(*,*) "hamsize:", hamsize
    call terminate
  end if

  allocate(isuppz(2*solsize))

  if(solsize == hamsize) then

    ! Get optimal work array lengths
    call workspacequery('a')
    allocate(work(lwork), rwork(lrwork), iwork(liwork))

    ! Diagonalize
    call zheevr('v', 'a', 'u', hamsize, ham, hamsize, vl, vu, il, iu, &
      & abstol, solsize, eval, evec, hamsize, isuppz, work, lwork, rwork,&
      & lrwork, iwork, liwork, info)

  else

    ! Smallest and largest eigenvalue indices
    il = 1
    iu = solsize

    ! Get optimal work array lengths
    call workspacequery('i')
    allocate(work(lwork), rwork(lrwork), iwork(liwork))

    ! Diagonalize
    call zheevr('v', 'i', 'u', hamsize, ham, hamsize, vl, vu, il, iu, &
      & abstol, solsize, eval, evec, hamsize, isuppz, work, lwork, rwork,&
      & lrwork, iwork, liwork, info)

  end if

  if(info .ne. 0) then
   write(*,*)
   write(*, '("Error(bsesoldiag): zheevx returned non-zero info:", i6)') info
   write(*,*)
   call terminate
  end if

  deallocate(isuppz)
  deallocate(work, rwork, iwork)

  contains

    subroutine workspacequery(rangetype)

      character(1), intent(in) :: rangetype

      lwork=-1
      lrwork=-1
      liwork=-1

      allocate(work(1), rwork(1), iwork(1))

      call zheevr('v', rangetype, 'u', hamsize, ham, hamsize, vl, vu, il, iu, &
        & abstol, solsize, eval, evec, hamsize, isuppz, work, lwork, rwork,&
        & lrwork, iwork, liwork, info)

      ! Adjust workspace
      lwork=int(work(1))
      lrwork=int(rwork(1))
      liwork=int(iwork(1))

      deallocate(work, rwork, iwork)

    end subroutine workspacequery

end subroutine bsesoldiag
!EOC

