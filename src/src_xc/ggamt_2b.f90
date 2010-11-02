

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_2b
! !INTERFACE:


subroutine ggamt_2b(is, g2rho, gvrho, vx, vc, dxdg2, dcdg2)
! !USES:
use modinput
use mod_Gvector
use mod_muffin_tin
use mod_SHT
use mod_atoms


! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_2b}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(in) :: g2rho(lmmaxvr, nrmtmax)
real(8), intent(in) :: gvrho(lmmaxvr, nrmtmax, 3)
real(8), intent(inout) :: vx(lmmaxvr, nrmtmax)
real(8), intent(inout) :: vc(lmmaxvr, nrmtmax)
real(8), intent(in) :: dxdg2(lmmaxvr, nrmtmax)
real(8), intent(in) :: dcdg2(lmmaxvr, nrmtmax)
! local variables
integer::nr, i
! allocatable arrays
real(8), allocatable :: rfmt1(:, :), rfmt2(:, :)
real(8), allocatable :: grfmt(:, :, :)
allocate(rfmt1(lmmaxvr, nrmtmax), rfmt2(lmmaxvr, nrmtmax))
allocate(grfmt(lmmaxvr, nrmtmax, 3))
nr=nrmt(is)
!------------------!
!     exchange     !
!------------------!
! convert dxdg2 to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dxdg2, lmmaxvr, 0.d0, &
 rfmt1, lmmaxvr)
! compute grad dxdg2
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dxdg2).(grad rho) in spherical coordinates
rfmt1(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt2, lmmaxvr)
  rfmt1(:, 1:nr)=rfmt1(:, 1:nr)+rfmt2(:, 1:nr)*gvrho(:, 1:nr, i)
end do
vx(:, 1:nr)=vx(:, 1:nr)-2.d0*(rfmt1(:, 1:nr)+dxdg2(:, 1:nr)*g2rho(:, 1:nr))
!---------------------!
!     correlation     !
!---------------------!
! convert dcdg2 to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dcdg2, lmmaxvr, 0.d0, &
 rfmt1, lmmaxvr)
! compute grad dcdg2
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dcdg2).(grad rho) in spherical coordinates
rfmt1(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt2, lmmaxvr)
  rfmt1(:, 1:nr)=rfmt1(:, 1:nr)+rfmt2(:, 1:nr)*gvrho(:, 1:nr, i)
end do
vc(:, 1:nr)=vc(:, 1:nr)-2.d0*(rfmt1(:, 1:nr)+dcdg2(:, 1:nr)*g2rho(:, 1:nr))
deallocate(rfmt1, rfmt2, grfmt)
return
end subroutine
!EOC
