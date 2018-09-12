

! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggair_2b
! !INTERFACE:


subroutine gv2xir(grho, vx, v2xsr)
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
real(8), intent(in) :: grho(ngrtot)
real(8), intent(inout) :: vx(ngrtot)
real(8), intent(in) :: v2xsr(ngrtot)
! local variables
integer::ig, ifg, i
real(8) :: g2rho(ngrtot)
real(8) :: gvrho(ngrtot, 3)
real(8) :: grho2(ngrtot)
! allocatable arrays
real(8), allocatable :: rfir(:)
complex(8), allocatable :: zfft1(:), zfft2(:)
allocate(rfir(ngrtot))
allocate(zfft1(ngrtot), zfft2(ngrtot))
call ggair_2a(g2rho,gvrho,grho2)
! compute grad dxdg2
zfft1(:)=v2xsr(:)
call zfftifc(3, ngrid, -1, zfft1)
! (grad dxdg2).(grad rho)
rfir(:)=0.d0
!!add gvrho
do i=1, 3
  zfft2(:)=0.d0
  do ig=1, ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i, ig)*cmplx(-aimag(zfft1(ifg)), dble(zfft1(ifg)), 8)
  end do
  call zfftifc(3, ngrid, 1, zfft2)
  rfir(:)=rfir(:)+dble(zfft2(:))*(gvrho(:, i))/grho(:)
end do
vx(:)=vx(:)-rfir(:)
return
end subroutine

