
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
integer isym,lspl,ilspl,lspn
integer i,md,sym(3,3),iv(3)
integer ig,ifg,jg,jfg
real(8) sc(3,3),vtc(3),t1
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
! inverse rotation required for rotation of G-vectors
  ilspl=isymlat(lspl)
  sym(:,:)=symlat(:,:,ilspl)
! global spin proper rotation in Cartesian coordinates
  lspn=lspnsymc(isym)
  md=symlatd(lspn)
  sc(:,:)=dble(md)*symlatc(:,:,lspn)
  do ig=1,ngv
    ifg=igfft(ig)
! multiply the transpose of the inverse symmetry matrix with the G-vector
    iv(1)=sym(1,1)*ivg(1,ig)+sym(2,1)*ivg(2,ig)+sym(3,1)*ivg(3,ig)
    iv(2)=sym(1,2)*ivg(1,ig)+sym(2,2)*ivg(2,ig)+sym(3,2)*ivg(3,ig)
    iv(3)=sym(1,3)*ivg(1,ig)+sym(2,3)*ivg(2,ig)+sym(3,3)*ivg(3,ig)
    iv(:)=modulo(iv(:)-intgv(:,1),ngrid(:))+intgv(:,1)
    jg=ivgig(iv(1),iv(2),iv(3))
    jfg=igfft(jg)
! complex phase factor for translation
    t1=-dot_product(vgc(:,ig),vtc(:))
    zt1=cmplx(cos(t1),sin(t1),8)
! translation, spatial rotation and global spin rotation
    if (lspn.eq.1) then
! global spin symmetry is the identity
      zfft2(jfg,:)=zfft2(jfg,:)+zt1*zfft1(ifg,:)
    else
      if (ndmag.eq.3) then
! non-collinear case
        zv(1)=sc(1,1)*zfft1(ifg,1)+sc(1,2)*zfft1(ifg,2)+sc(1,3)*zfft1(ifg,3)
        zv(2)=sc(2,1)*zfft1(ifg,1)+sc(2,2)*zfft1(ifg,2)+sc(2,3)*zfft1(ifg,3)
        zv(3)=sc(3,1)*zfft1(ifg,1)+sc(3,2)*zfft1(ifg,2)+sc(3,3)*zfft1(ifg,3)
        zfft2(jfg,:)=zfft2(jfg,:)+zt1*zv(:)
      else
! collinear case
        zfft2(jfg,1)=zfft2(jfg,1)+sc(3,3)*zt1*zfft1(ifg,1)
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

