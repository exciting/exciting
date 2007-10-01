
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrf
! !INTERFACE:
subroutine symrf(lrstp,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rfmt  : real muffin-tin function (inout,real(lmmaxvr,nrmtmax,natmtot))
!   rfir  : real intersitial function (inout,real(ngrtot))
! !DESCRIPTION:
!   Symmetrises a real scalar function defined over the entire unit cell using
!   the full set of crystal symmetries. In the muffin-tin of a particular atom
!   the spherical harmonic coefficients of every equivlent atom are rotated and
!   averaged. The interstitial part of the function is first Fourier transformed
!   to $G$-space, and then averaged over each symmetry by rotating the Fourier
!   coefficients and multiplying them by a phase factor corresponding to the
!   symmetry translation.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(inout) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(inout) :: rfir(ngrtot)
! local variables
integer is,ia,ja,ias,jas
integer isym,lspl,sym(3,3)
integer ig,jg,ifg,jfg
integer iv(3),ir
real(8) vtc(3),t1,t2
complex(8) zt1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rfmt1(:,:,:),rfmt2(:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
allocate(rfmt1(lmmaxvr,nrmtmax,natmmax))
allocate(rfmt2(lmmaxvr,nrmtmax))
allocate(zfft1(ngrtot),zfft2(ngrtot))
t1=1.d0/dble(nsymcrys)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
! make a copy of the input function
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is),lrstp
      rfmt1(:,ir,ia)=rfmt(:,ir,ias)
    end do
  end do
  done(:)=.false.
! loop over atoms
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        rfmt(:,ir,ias)=0.d0
      end do
! loop over crystal symmetries
      do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
        lspl=lsplsymc(isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
        ja=ieqatom(ia,is,isym)
! apply the rotation to the muffin-tin function
        call symrfmt(lrstp,is,symlat(1,1,lspl),rfmt1(1,1,ja),rfmt2)
! accumulate in original function array
        do ir=1,nrmt(is),lrstp
          rfmt(:,ir,ias)=rfmt(:,ir,ias)+rfmt2(:,ir)
        end do
      end do
! normalise
      do ir=1,nrmt(is),lrstp
        rfmt(:,ir,ias)=t1*rfmt(:,ir,ias)
      end do
      done(ia)=.true.
! rotate into equivalent atoms
      do isym=1,nsymcrys
        ja=ieqatom(ia,is,isym)
        if (.not.done(ja)) then
          jas=idxas(ja,is)
          lspl=lsplsymc(isym)
! find inverse symmetry (which rotates atom ia into atom ja)
          call i3minv(symlat(1,1,lspl),sym)
! rotate symmetrised function into equivalent muffin-tin
          call symrfmt(lrstp,is,sym,rfmt(1,1,ias),rfmt(1,1,jas))
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
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
  do ig=1,ngvec
    ifg=igfft(ig)
    t2=-dot_product(vgc(:,ig),vtc(:))
! complex phase factor for translation
    zt1=cmplx(cos(t2),sin(t2),8)
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
rfir(:)=t1*dble(zfft2(:))
deallocate(rfmt1,rfmt2,zfft1,zfft2)
return
end subroutine
!EOC

