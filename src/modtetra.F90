
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtetra
  ! Variable names taken from the GW implementation into the EXCITING code
  ! version 0.9.52 by R. Gomez-Abal.
  implicit none

  !----------------------------!
  !     ordering variables     !
  !----------------------------!
  ! map from library k-point index to application k-point index
  integer, allocatable :: iktet2ik(:)
  ! reverse map
  integer, allocatable :: ik2iktet(:)

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

  subroutine geniktetmap(eps,nppt,ngridp,vploff,vpllib,vpl,ipmap)
    implicit none
    ! arguments
    real(8), intent(in) :: eps
    integer, intent(in) :: nppt
    integer, intent(in) :: ngridp(3)
    real(8), intent(in) :: vploff(3)
    real(8), intent(in) :: vpl(3,nppt)
    real(8), intent(in) :: vpllib(3,nppt)
    integer, intent(in) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
    ! local variables
    integer :: ip,ipd,iv(3)
    if (allocated(iktet2ik)) deallocate(iktet2ik)
    allocate(iktet2ik(nppt))
    if (allocated(ik2iktet)) deallocate(ik2iktet)
    allocate(ik2iktet(nppt))
    do ip=1,nppt
       ! grid coordinates of library k-point
       iv(:)=nint(vpllib(:,ip)*ngridp-vploff(:))
       ! index in default p-point set
       ipd=ipmap(iv(1),iv(2),iv(3))
       ! map from library to default
       iktet2ik(ip)=ipd
       ! reverse map
       ik2iktet(ipd)=ip
       ! check maps
!!$       if (sum(abs(vpl(:,ipd)-vpllib(:,ip))).ge.eps) then
!!$          write(*,*)
!!$          write(*,'("Error(modtetra:geniktetmap): k-point mapping between")')
!!$          write(*,'(" set of library and default set failed:")')
          write(*,'(" library k-point       :",3g18.10)') vpllib(:,ip)
          write(*,'(" mapped default k-point:",3g18.10)') vpl(:,ipd)
          write(*,'(" map from library to default:",2i8)') ip,ipd
          write(*,*)
!!$          stop
!!$       end if
    end do
  end subroutine geniktetmap

end module modtetra
