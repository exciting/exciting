
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvfir
subroutine symrvfir(ngv,rvfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngv   : number of G-vectors to be used for the Fourier space rotation
!           (in,integer)
!   rvfir : real interstitial vector function (inout,real(ngrtot,ndmag))
! !DESCRIPTION:
!   Symmetrises a real interstitial vector function. See routines {\tt symrvf}
!   and {\tt symrfir} for details.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngv
real(8), intent(inout) :: rvfir(ngrtot,ndmag)
! local variables
integer i,isym,lspl,lspn
integer ig,ifg,jg,jfg,iv(3)
real(8) s(3,3),vtc(3),t1
complex(8) zv(3),zt1
! allocatable arrays
complex(8), allocatable :: zfft1(:,:),zfft2(:,:)
allocate(zfft1(ngrtot,ndmag),zfft2(ngrtot,ndmag))
! Fourier transform vector function to G-space
do i=1,ndmag
  zfft1(:,i)=rvfir(:,i)
  call zfftifc(3,ngrid,-1,zfft1(1,i))
end do
zfft2(:,:)=0.d0
do isym=1,nsymcrys
! translation vector in Cartesian coordinates
  call r3mv(avec,vtlsymc(1,isym),vtc)
! index to spatial rotation lattice symmetry
  lspl=lsplsymc(isym)
! global spin rotation in Cartesian coordinates
  lspn=lspnsymc(isym)
  s(:,:)=dble(symlat(:,:,lspn))
  call r3mm(s,ainv,s)
  call r3mm(avec,s,s)
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
! translation, spatial rotation and global spin rotation
    if (lspn.eq.1) then
! global spin symmetry is the identity
      zfft2(jfg,:)=zfft2(jfg,:)+zt1*zfft1(ifg,:)
    else
      if (ndmag.eq.3) then
! non-collinear case
        zv(1)=s(1,1)*zfft1(ifg,1)+s(1,2)*zfft1(ifg,2)+s(1,3)*zfft1(ifg,3)
        zv(2)=s(2,1)*zfft1(ifg,1)+s(2,2)*zfft1(ifg,2)+s(2,3)*zfft1(ifg,3)
        zv(3)=s(3,1)*zfft1(ifg,1)+s(3,2)*zfft1(ifg,2)+s(3,3)*zfft1(ifg,3)
        zfft2(jfg,:)=zfft2(jfg,:)+zt1*zv(:)
      else
! collinear case
        zfft2(jfg,1)=zfft2(jfg,1)+s(3,3)*zt1*zfft1(ifg,1)
      end if
    end if
  end do
end do
! Fourier transform to real-space and normalise
t1=1.d0/dble(nsymcrys)
do i=1,ndmag
  call zfftifc(3,ngrid,1,zfft2(1,i))
  rvfir(:,i)=t1*dble(zfft2(:,i))
end do
deallocate(zfft1,zfft2)
return
end subroutine
!EOC

