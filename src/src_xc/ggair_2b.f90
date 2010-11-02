

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2b
! !INTERFACE:


subroutine ggair_2b(g2rho, gvrho, vx, vc, dxdg2, dcdg2)
! !USES:
use modinput
use mod_Gvector
use mod_potential_and_density
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_2b}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: g2rho(ngrtot)
real(8), intent(in) :: gvrho(ngrtot, 3)
real(8), intent(inout) :: vx(ngrtot)
real(8), intent(inout) :: vc(ngrtot)
real(8), intent(in) :: dxdg2(ngrtot)
real(8), intent(in) :: dcdg2(ngrtot)
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
! compute grad dxdg2
zfft1(:)=dxdg2(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dxdg2).(grad rho)
rfir(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvrho(:, i)
end do
vx(:)=vx(:)-2.d0*(rfir(:)+dxdg2(:)*g2rho(:))
!---------------------!
!     correlation     !
!---------------------!
! compute grad dcdg2
zfft1(:)=dcdg2(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dcdg2).(grad rho)
rfir(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*gvrho(:, i)
end do
vc(:)=vc(:)-2.d0*(rfir(:)+dcdg2(:)*g2rho(:))
deallocate(rfir, zfft1, zfft2)
return
end subroutine
!EOC
