
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtetra
  ! Variable names taken from the GW implementation into the EXCITING code
  ! version 0.9.52 by R. Gomez-Abal.
  implicit none

  !--------------------------------------!
  !     tetrahedron method variables     !
  !--------------------------------------!
  ! true if tetrahedron method is used
  logical :: tetra
  ! integer k-point offset
  integer(4) :: ikloff(3)
  ! k-point offset divisor
  integer(4) :: dkloff
  ! k-points common divisor
  integer(4) :: dvk
  ! Number of tetrahedra
  integer(4) :: ntet
  ! index of the k-points corresponding to the nodes of each tetrahedron
  integer(4), allocatable :: tnodes(:,:)
  ! weight of each tetrahedron.
  integer(4), allocatable :: wtet(:)
  ! volume of the tetrahedra relative to the BZ volume
  real(8) :: tvol
  ! parameter specifying smalles diagonal in generic tetrahedron
  integer(4) :: mnd

  !---------------------------------!
  !     q-dependent convolution     !
  !---------------------------------!
  ! number of the tetrahedra linked to by the corresponding q vector
  integer(4), allocatable :: link(:), kqid(:,:)
  ! q-points common divisor
  integer(4) dvq

contains

!BOP
! !ROUTINE: rtorat
! !INTERFACE:
  subroutine rtorat(eps,n,x,k,div)
! !DESCRIPTION: 
!   This subroutine factorizes the real coordinates of a vector {\bf x}.
!   The output is an integer vector {\bf k}, such that $k(i)/{\rm div}=x(i)$
!   and
!   $$ |x(i)-k(i)/{\rm div}| < {\rm eps} $$.
!
! !REVISION HISTORY:
!   Created July 2008 by Sagmeister
!EOP
!BOC
    implicit none
    ! arguments
    real(8), intent(in) :: eps
    integer(4), intent(in) :: n
    real(8), intent(in) :: x(n)
    integer(4), intent(out) :: div
    integer(4), intent(out) :: k(n)
    ! local variables
    integer :: maxint
    real(8) :: dx
    maxint=nint(1.d0/eps)
    do div=1,maxint
       k(:)=nint(dble(div)*x(:))
       dx=maxval(abs(dble(k)/dble(div)-x))
       if (dx.lt.eps) exit
    end do
    if (dx.ge.eps) then
       write(*,*)
       write(*,'("Error(modtetra:rtorat): factorization failed")')
       write(*,'(" maximum integer :",i12)') maxint
       write(*,'(" tolerance       :",g18.10)') eps
       write(*,'(" deviation       :",g18.10)') dx
       write(*,*)
       stop
    end if
    if (dx.gt.1.d-12) then
       write(*,*)
       write(*,'("Warning(modtetra:rtorat): small deviation in factorization")')
       write(*,'(" maximum deviation :",g18.10)') dx
       write(*,*)
    end if
  end subroutine rtorat
!EOC

  subroutine r3fraction(r,n,d)
    implicit none
    ! arguments
    real(8), intent(in) :: r(3)
    integer, intent(out) :: n(3),d
    ! parameters
    real(8), parameter :: eps=1.d-5,eps2=1.d-3
    call rtorat(eps,3,r,n,d)
    ! check factorization
    if ((sum(abs(r)).lt.eps2).and.(sum(abs(r)).gt.0.d0)) then
       write(*,*)
       write(*,'("Warning(modtetra:r3fraction): very small offset:")')
       write(*,'(" kgen and related routines might fail")')
       write(*,*)
    end if
  end subroutine r3fraction

end module modtetra
