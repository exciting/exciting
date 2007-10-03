
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrfir
subroutine symrfir(ngv,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngv  : number of G-vectors to be used for the Fourier space rotation
!          (in,integer)
!   rfir : real intersitial function (inout,real(ngrtot))
! !DESCRIPTION:
!   Symmetrises a real scalar interstitial function. The function is first
!   Fourier transformed to $G$-space, and then averaged over each symmetry by
!   rotating the Fourier coefficients and multiplying them by a phase factor
!   corresponding to the symmetry translation.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngv
real(8), intent(inout) :: rfir(ngrtot)
! local variables
integer isym,lspl,iv(3)
integer ig,jg,ifg,jfg
real(8) vtc(3),t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(zfft1(ngrtot),zfft2(ngrtot))
! Fourier transform function to G-space
zfft1(:)=rfir(:)
call zfftifc(3,ngrid,-1,zfft1)
zfft2(:)=0.d0
! loop over crystal symmetries
do isym=1,nsymcrys
! translation in Cartesian coordinates
  call r3mv(avec,vtlsymc(1,isym),vtc)
! index to lattice symmetry of spatial rotation
  lspl=lsplsymc(isym)
  do ig=1,ngv
    ifg=igfft(ig)
    t1=-dot_product(vgc(:,ig),vtc(:))
! complex phase factor for translation
    zt1=cmplx(cos(t1),sin(t1),8)
! multiply the transpose of the symmetry matrix with the G-vector
    iv(1)=symlat(1,1,lspl)*ivg(1,ig) &
         +symlat(2,1,lspl)*ivg(2,ig) &
         +symlat(3,1,lspl)*ivg(3,ig)
    iv(2)=symlat(1,2,lspl)*ivg(1,ig) &
         +symlat(2,2,lspl)*ivg(2,ig) &
         +symlat(3,2,lspl)*ivg(3,ig)
    iv(3)=symlat(1,3,lspl)*ivg(1,ig) &
         +symlat(2,3,lspl)*ivg(2,ig) &
         +symlat(3,3,lspl)*ivg(3,ig)
    iv(:)=modulo(iv(:)-intgv(:,1),ngrid(:))+intgv(:,1)
    jg=ivgig(iv(1),iv(2),iv(3))
    jfg=igfft(jg)
    zfft2(jfg)=zfft2(jfg)+zt1*zfft1(ifg)
  end do
end do
! Fourier transform to real-space and normalise
call zfftifc(3,ngrid,1,zfft2)
t1=1.d0/dble(nsymcrys)
rfir(:)=t1*dble(zfft2(:))
deallocate(zfft1,zfft2)
return
end subroutine
!EOC

