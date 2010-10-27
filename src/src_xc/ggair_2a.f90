

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2a
! !INTERFACE:


subroutine ggair_2a(g2rho, gvrho, grho2)
! !USES:
use modinput
use mod_Gvector
use mod_potential_and_density
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggair\_sp\_2a}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: g2rho(ngrtot)
real(8), intent(out) :: gvrho(ngrtot, 3)
real(8), intent(out) :: grho2(ngrtot)
! local variables
integer::i, ig, ifg
! allocatable arrays
complex(8), allocatable :: zfft1(:), zfft2(:)
allocate(zfft1(ngrtot), zfft2(ngrtot))
! Fourier transform density to G-space
zfft1(:)=rhoir(:)
call zfftifc(3, ngrid, -1, zfft1)
! compute grad^2 rho
zfft2(:)=0.d0
do ig=1, ngvec
  ifg=igfft(ig)
  zfft2(ifg)=-(gc(ig)**2)*zfft1(ifg)
end do
call zfftifc(3, ngrid, 1, zfft2)
g2rho(:)=dble(zfft2(:))
! compute grad rho and (grad rho)^2
grho2(:)=0.d0
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  gvrho(:, i)=dble(zfft2(:))
  grho2(:)=grho2(:)+gvrho(:, i)**2
end do
deallocate(zfft1, zfft2)
return
end subroutine
!EOC
