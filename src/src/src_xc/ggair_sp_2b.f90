

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_sp_2b
! !INTERFACE:


subroutine ggair_sp_2b(g2up, g2dn, gvup, gvdn, vxup, vxdn, vcup, vcdn, dxdgu2, dxdgd2, &
 dxdgud, dcdgu2, dcdgd2, dcdgud)
use modinput
! !USES:
use mod_Gvector
! !DESCRIPTION:
!   Post processing step of interstitial gradients for GGA type 2. See routine
!   {\tt ggamt\_sp\_2a} for full details.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
real(8), intent(in) :: g2up(ngrtot)
real(8), intent(in) :: g2dn(ngrtot)
real(8), intent(in) :: gvup(ngrtot, 3)
real(8), intent(in) :: gvdn(ngrtot, 3)
real(8), intent(inout) :: vxup(ngrtot)
real(8), intent(inout) :: vxdn(ngrtot)
real(8), intent(inout) :: vcup(ngrtot)
real(8), intent(inout) :: vcdn(ngrtot)
real(8), intent(in) :: dxdgu2(ngrtot)
real(8), intent(in) :: dxdgd2(ngrtot)
real(8), intent(in) :: dxdgud(ngrtot)
real(8), intent(in) :: dcdgu2(ngrtot)
real(8), intent(in) :: dcdgd2(ngrtot)
real(8), intent(in) :: dcdgud(ngrtot)
! local variables
integer::ig, ifg, i
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfft1(:), zfft2(:)
allocate(rfir(ngrtot))
allocate(zfft1(ngrtot), zfft2(ngrtot))
!------------------!
!     exchange     !
!------------------!
! compute grad dxdgu2
zfft1(:)=dxdgu2(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dxdgu2).(grad rhoup)
rfir(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvup(:, i)
end do
vxup(:)=vxup(:)-2.d0*(rfir(:)+dxdgu2(:)*g2up(:))-dxdgud(:)*g2dn(:)
! compute grad dxdgd2
zfft1(:)=dxdgd2(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dxdgd2).(grad rhodn)
rfir(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvdn(:, i)
end do
vxdn(:)=vxdn(:)-2.d0*(rfir(:)+dxdgd2(:)*g2dn(:))-dxdgud(:)*g2up(:)
! compute grad dxdgud
zfft1(:)=dxdgud(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dxdgud).(grad rhodn) and (grad dxdgud).(grad rhoup)
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  vxup(:)=vxup(:)-dble(zfft2(:))*gvdn(:, i)
  vxdn(:)=vxdn(:)-dble(zfft2(:))*gvup(:, i)
end do
!---------------------!
!     correlation     !
!---------------------!
! compute grad dcdgu2
zfft1(:)=dcdgu2(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dcdgu2).(grad rhoup)
rfir(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvup(:, i)
end do
vcup(:)=vcup(:)-2.d0*(rfir(:)+dcdgu2(:)*g2up(:))-dcdgud(:)*g2dn(:)
! compute grad dcdgd2
zfft1(:)=dcdgd2(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dcdgd2).(grad rhodn)
rfir(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvdn(:, i)
end do
vcdn(:)=vcdn(:)-2.d0*(rfir(:)+dcdgd2(:)*g2dn(:))-dcdgud(:)*g2up(:)
! compute grad dcdgud
zfft1(:)=dcdgud(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dcdgud).(grad rhodn) and (grad dcdgud).(grad rhoup)
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  vcup(:)=vcup(:)-dble(zfft2(:))*gvdn(:, i)
  vcdn(:)=vcdn(:)-dble(zfft2(:))*gvup(:, i)
end do
deallocate(rfir, zfft1, zfft2)
return
end subroutine
!EOC
