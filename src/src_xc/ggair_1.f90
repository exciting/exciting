

! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_1
! !INTERFACE:


subroutine ggair_1(grho, g2rho, g3rho)
! !USES:
use modinput
use mod_Gvector
use mod_potential_and_density
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_1}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: grho(ngrtot)
real(8), intent(out) :: g2rho(ngrtot)
real(8), intent(out) :: g3rho(ngrtot)
! local variables
integer::i, ig, ifg
! allocatable arrays
real(8), allocatable :: gvrho(:, :)
complex(8), allocatable :: zfft1(:), zfft2(:)
allocate(gvrho(ngrtot, 3))
allocate(zfft1(ngrtot), zfft2(ngrtot))
zfft1(:)=rhoir(:)
call zfftifc(3, ngrid, -1, zfft1)
! |grad rho|
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  gvrho(:, i)=dble(zfft2(:))
end do
grho(:)=sqrt(gvrho(:, 1)**2+gvrho(:, 2)**2+gvrho(:, 3)**2)
! grad^2 rho
zfft2(:)=0.d0
do ig=1, ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3, ngrid, 1, zfft2)
g2rho(:)=dble(zfft2(:))
! (grad rho).(grad |grad rho|)
zfft1(:)=grho(:)
call zfftifc(3, ngrid, -1, zfft1)
g3rho(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  g3rho(:)=g3rho(:)+gvrho(:, i)*dble(zfft2(:))
end do
deallocate(gvrho, zfft1, zfft2)
return
end subroutine
!EOC
