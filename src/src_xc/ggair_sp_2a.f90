

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_sp_2a
! !INTERFACE:


subroutine ggair_sp_2a(rhoup, rhodn, g2up, g2dn, gvup, gvdn, gup2, gdn2, gupdn)
! !USES:
use modinput
use mod_Gvector
! !DESCRIPTION:
!   Computes the interstitial gradients $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho^{\uparrow}$,
!   $\nabla\rho^{\downarrow}$, $(\nabla\rho^{\uparrow})^2$,
!   $(\nabla\rho^{\downarrow})^2$ and
!   $\nabla\rho^{\uparrow}\cdot\nabla\rho^{\downarrow}$. These are used for GGA
!   functionals of type 2. See {\tt ggamt\_sp\_2a} for full details.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rhoup(ngrtot)
real(8), intent(in) :: rhodn(ngrtot)
real(8), intent(out) :: g2up(ngrtot)
real(8), intent(out) :: g2dn(ngrtot)
real(8), intent(out) :: gvup(ngrtot, 3)
real(8), intent(out) :: gvdn(ngrtot, 3)
real(8), intent(out) :: gup2(ngrtot)
real(8), intent(out) :: gdn2(ngrtot)
real(8), intent(out) :: gupdn(ngrtot)
! local variables
integer::ig, ifg, i
! allocatable arrays
complex(8), allocatable :: zfft1(:)
complex(8), allocatable :: zfft2(:)
allocate(zfft1(ngrtot))
allocate(zfft2(ngrtot))
!----------------!
!     rho up     !
!----------------!
zfft1(:)=rhoup(:)
call zfftifc(3, ngrid, -1, zfft1)
! compute grad^2 rhoup
zfft2(:)=0.d0
do ig=1, ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3, ngrid, 1, zfft2)
g2up(:)=dble(zfft2(:))
! grad rhoup
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  gvup(:, i)=dble(zfft2(:))
end do
! (grad rhoup)^2
gup2(:)=gvup(:, 1)**2+gvup(:, 2)**2+gvup(:, 3)**2
!------------------!
!     rho down     !
!------------------!
zfft1(:)=rhodn(:)
call zfftifc(3, ngrid, -1, zfft1)
! compute grad^2 rhodn
zfft2(:)=0.d0
do ig=1, ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3, ngrid, 1, zfft2)
g2dn(:)=dble(zfft2(:))
! grad rhodn
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  gvdn(:, i)=dble(zfft2(:))
end do
! (grad rhodn)^2
gdn2(:)=gvdn(:, 1)**2+gvdn(:, 2)**2+gvdn(:, 3)**2
! (grad rhoup).(grad rhodn)
gupdn(:) = gvup(:, 1) * gvdn(:, 1) &
	+gvup(:, 2) * gvdn(:, 2) &
	+gvup(:, 3) * gvdn(:, 3)
deallocate(zfft1, zfft2)
return
end subroutine
!EOC
