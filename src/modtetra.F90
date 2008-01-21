
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtetra
  ! variable names taken from the GW implementation into the EXCITING code
  ! version 0.9.52 done by R. Gomez-Abal
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

  subroutine r3fraction(r,n,d)
    implicit none
    ! arguments
    real(8), intent(in) :: r(3)
    integer, intent(out) :: n(3),d
    ! parameters
    real(8), parameter :: epst=1.d-6
    ! call to libbzint-routine
    call factorize(3,r,n,d)
    ! check factorization
    if (any(abs(dble(n)/dble(d)-r).gt.epst)) then
       write(*,*)
       write(*,'("Error(modtetra:r3fraction): factorization failed:")')
       write(*,'(" vector                   :",3g18.10)') r
       write(*,'(" vector from factorization:",3g18.10)') dble(n)/dble(d)
       write(*,*)
       stop
    end if
  end subroutine r3fraction

end module modtetra
