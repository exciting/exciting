

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_sp_2b
! !INTERFACE:


subroutine ggamt_sp_2b(is, g2up, g2dn, gvup, gvdn, vxup, vxdn, vcup, vcdn, dxdgu2, &
 dxdgd2, dxdgud, dcdgu2, dcdgd2, dcdgud)

! !USES:
use modinput
use mod_Gvector
use mod_muffin_tin
use mod_SHT
use mod_atoms
! !DESCRIPTION:
!   Post processing step of muffin-tin gradients for GGA type 2. See routine
!   {\tt ggamt\_sp\_2a} for full details.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
 ! arguments
integer, intent(in) :: is
real(8), intent(in) :: g2up(lmmaxvr, nrmtmax)
real(8), intent(in) :: g2dn(lmmaxvr, nrmtmax)
real(8), intent(in) :: gvup(lmmaxvr, nrmtmax, 3)
real(8), intent(in) :: gvdn(lmmaxvr, nrmtmax, 3)
real(8), intent(inout) :: vxup(lmmaxvr, nrmtmax)
real(8), intent(inout) :: vxdn(lmmaxvr, nrmtmax)
real(8), intent(inout) :: vcup(lmmaxvr, nrmtmax)
real(8), intent(inout) :: vcdn(lmmaxvr, nrmtmax)
real(8), intent(in) :: dxdgu2(lmmaxvr, nrmtmax)
real(8), intent(in) :: dxdgd2(lmmaxvr, nrmtmax)
real(8), intent(in) :: dxdgud(lmmaxvr, nrmtmax)
real(8), intent(in) :: dcdgu2(lmmaxvr, nrmtmax)
real(8), intent(in) :: dcdgd2(lmmaxvr, nrmtmax)
real(8), intent(in) :: dcdgud(lmmaxvr, nrmtmax)
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
! convert dxdgu2 to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dxdgu2, lmmaxvr, &
 0.d0, rfmt1, lmmaxvr)
! compute grad dxdgu2
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dxdgu2).(grad rhoup) in spherical coordinates
rfmt1(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt2, lmmaxvr)
  rfmt1(:, 1:nr)=rfmt1(:, 1:nr)+rfmt2(:, 1:nr)*gvup(:, 1:nr, i)
end do
vxup(:, 1:nr) = vxup(:, 1:nr) - 2.d0 * (rfmt1(:, 1:nr) + dxdgu2(:, 1:nr) * g2up(:, 1:nr)) &
 -dxdgud(:, 1:nr) * g2dn(:, 1:nr)
! convert dxdgd2 to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dxdgd2, lmmaxvr, &
 0.d0, rfmt1, lmmaxvr)
! compute grad dxdgd2
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dxdgd2).(grad rhodn) in spherical coordinates
rfmt1(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt2, lmmaxvr)
  rfmt1(:, 1:nr)=rfmt1(:, 1:nr)+rfmt2(:, 1:nr)*gvdn(:, 1:nr, i)
end do
vxdn(:, 1:nr) = vxdn(:, 1:nr) - 2.d0 * (rfmt1(:, 1:nr) + dxdgd2(:, 1:nr) * g2dn(:, 1:nr)) &
 -dxdgud(:, 1:nr) * g2up(:, 1:nr)
! convert dxdgud to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dxdgud, lmmaxvr, &
 0.d0, rfmt1, lmmaxvr)
! compute grad dxdgud
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dxdgud).(grad rhodn) and (grad dxdgud).(grad rhoup)
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt1, lmmaxvr)
  vxup(:, 1:nr)=vxup(:, 1:nr)-rfmt1(:, 1:nr)*gvdn(:, 1:nr, i)
  vxdn(:, 1:nr)=vxdn(:, 1:nr)-rfmt1(:, 1:nr)*gvup(:, 1:nr, i)
end do
!---------------------!
!     correlation     !
!---------------------!
! convert dcdgu2 to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dcdgu2, lmmaxvr, &
 0.d0, rfmt1, lmmaxvr)
! compute grad dcdgu2
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dcdgu2).(grad rhoup) in spherical coordinates
rfmt1(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt2, lmmaxvr)
  rfmt1(:, 1:nr)=rfmt1(:, 1:nr)+rfmt2(:, 1:nr)*gvup(:, 1:nr, i)
end do
vcup(:, 1:nr) = vcup(:, 1:nr) - 2.d0 * (rfmt1(:, 1:nr) + dcdgu2(:, 1:nr) * g2up(:, 1:nr)) &
 -dcdgud(:, 1:nr) * g2dn(:, 1:nr)
! convert dcdgd2 to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dcdgd2, lmmaxvr, &
 0.d0, rfmt1, lmmaxvr)
! compute grad dcdgd2
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dcdgd2).(grad rhodn) in spherical coordinates
rfmt1(:, 1:nr)=0.d0
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt2, lmmaxvr)
  rfmt1(:, 1:nr)=rfmt1(:, 1:nr)+rfmt2(:, 1:nr)*gvdn(:, 1:nr, i)
end do
vcdn(:, 1:nr) = vcdn(:, 1:nr) - 2.d0 * (rfmt1(:, 1:nr) + dcdgd2(:, 1:nr) * g2dn(:, 1:nr)) &
 -dcdgud(:, 1:nr) * g2up(:, 1:nr)
! convert dcdgud to spherical harmonics
call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rfshtvr, lmmaxvr, dcdgud, lmmaxvr, &
 0.d0, rfmt1, lmmaxvr)
! compute grad dcdgud
call gradrfmt(input%groundstate%lmaxvr, nr, spr(:, is), lmmaxvr, nrmtmax, rfmt1, grfmt)
! (grad dcdgud).(grad rhodn) and (grad dcdgud).(grad rhoup)
do i=1, 3
  call dgemm('N', 'N', lmmaxvr, nr, lmmaxvr, 1.d0, rbshtvr, lmmaxvr, grfmt(:, :, i), &
   lmmaxvr, 0.d0, rfmt1, lmmaxvr)
  vcup(:, 1:nr)=vcup(:, 1:nr)-rfmt1(:, 1:nr)*gvdn(:, 1:nr, i)
  vcdn(:, 1:nr)=vcdn(:, 1:nr)-rfmt1(:, 1:nr)*gvup(:, 1:nr, i)
end do
deallocate(rfmt1, rfmt2, grfmt)
return
end subroutine
!EOC
