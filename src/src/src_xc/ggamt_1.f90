

! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_1
! !INTERFACE:


subroutine ggamt_1(is, ia, grho, g2rho, g3rho)
! !USES:
use modinput
use mod_Gvector
use mod_muffin_tin
use mod_SHT
use mod_atoms
use mod_potential_and_density
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_1}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD)
!   Modified Mai 2012 (UW)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: ia
real(8), intent(out) :: grho(lmmaxvr, nrmtmax)
real(8), intent(out) :: g2rho(lmmaxvr, nrmtmax)
real(8), intent(out) :: g3rho(lmmaxvr, nrmtmax)
! local variables
integer::ias, nr, i
! allocatable arrays
real(8), allocatable :: grfmt(:, :, :), gvrho(:, :, :), rfmt(:, :)
allocate(grfmt(lmmaxvr, nrmtmax, 3), gvrho(lmmaxvr, nrmtmax, 3))
allocate(rfmt(lmmaxvr, nrmtmax))
ias=idxas(ia, is)
nr=nrmt(is)
! |grad rho|
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr,&
& nrmtmax, rhomt(:, :, ias), grfmt)
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, gvrho(:, :, i), lmmaxvr)
end do
grho(:, 1:nr)=sqrt(gvrho(:, 1:nr, 1)**2+gvrho(:, 1:nr, 2)**2+gvrho(:, 1:nr, 3)**2)
! grad^2 rho
call grad2rfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, rhomt(:, :, ias), rfmt)
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, rfmt, lmmaxvr, 0.d0, &
 g2rho, lmmaxvr)
! (grad rho).(grad |grad rho|)
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr,&
& grho, lmmaxvr, 0.d0, &
 rfmt, lmmaxvr)
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt, grfmt)
g3rho(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt, lmmaxvr)
  g3rho(:, 1:nr)=g3rho(:, 1:nr)+gvrho(:, 1:nr, i)*rfmt(:, 1:nr)
end do
deallocate(grfmt, gvrho, rfmt)
return
end subroutine
!EOC
