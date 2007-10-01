
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrvf(lrstp,rvfmt,rvfir)
use modmain
implicit none
! arguments
integer, intent(in) :: lrstp
real(8), intent(inout) :: rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag)
real(8), intent(inout) :: rvfir(ngrtot,ndmag)
! local variables
integer is,ia,ja,ias,jas
integer isym,lspl,lspn,sym(3,3)
integer ig,jg,ifg,jfg
integer iv(3),ir,lm,i
real(8) s(3,3),vtc(3),v(3),t1,t2
complex(8) zv(3),zt1
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rvfmt1(:,:,:,:),rvfmt2(:,:,:)
complex(8), allocatable :: zfft1(:,:),zfft2(:,:)
allocate(rvfmt1(lmmaxvr,nrmtmax,natmmax,ndmag))
allocate(rvfmt2(lmmaxvr,nrmtmax,ndmag))
allocate(zfft1(ngrtot,ndmag),zfft2(ngrtot,ndmag))
t1=1.d0/dble(nsymcrys)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
! make copy of vector field for all atoms of current species
  do i=1,ndmag
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        rvfmt1(:,ir,ia,i)=rvfmt(:,ir,ias,i)
      end do
    end do
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        rvfmt(:,ir,ias,:)=0.d0
      end do
! begin loop over crystal symmetries
      do isym=1,nsymcrys
! equivalent atom
        ja=ieqatom(ia,is,isym)
! parallel transport of vector field
        lspl=lsplsymc(isym)
        do i=1,ndmag
          call symrfmt(lrstp,is,symlat(1,1,lspl),rvfmt1(1,1,ja,i),rvfmt2(1,1,i))
        end do
! global spin rotation matrix in Cartesian coordinates
        lspn=lspnsymc(isym)
        s(:,:)=dble(symlat(:,:,lspn))
        call r3mm(s,ainv,s)
        call r3mm(avec,s,s)
! global spin rotation of vector field
        if (ndmag.eq.3) then
! non-collinear case
          do ir=1,nrmt(is),lrstp
            do lm=1,lmmaxvr
              v(1)=s(1,1)*rvfmt2(lm,ir,1) &
                  +s(1,2)*rvfmt2(lm,ir,2) &
                  +s(1,3)*rvfmt2(lm,ir,3)
              v(2)=s(2,1)*rvfmt2(lm,ir,1) &
                  +s(2,2)*rvfmt2(lm,ir,2) &
                  +s(2,3)*rvfmt2(lm,ir,3)
              v(3)=s(3,1)*rvfmt2(lm,ir,1) &
                  +s(3,2)*rvfmt2(lm,ir,2) &
                  +s(3,3)*rvfmt2(lm,ir,3)
              rvfmt(lm,ir,ias,:)=rvfmt(lm,ir,ias,:)+v(:)
            end do
          end do
        else
! collinear case
          do ir=1,nrmt(is),lrstp
            do lm=1,lmmaxvr
              rvfmt(lm,ir,ias,1)=rvfmt(lm,ir,ias,1)+s(3,3)*rvfmt2(lm,ir,1)
            end do
          end do
        end if
! end loop over crystal symmetries
      end do
! normalise
      do ir=1,nrmt(is),lrstp
        rvfmt(:,ir,ias,:)=t1*rvfmt(:,ir,ias,:)
      end do
! mark atom as done
      done(ia)=.true.
! rotate into equivalent atoms
      do isym=1,nsymcrys
        ja=ieqatom(ia,is,isym)
        if (.not.done(ja)) then
          jas=idxas(ja,is)
! parallel transport of vector field (using operation inverse)
          lspl=lsplsymc(isym)
          call i3minv(symlat(1,1,lspl),sym)
          do i=1,ndmag
            call symrfmt(lrstp,is,sym,rvfmt(1,1,ias,i),rvfmt(1,1,jas,i))
          end do
! inverse of global rotation matrix in Cartesian coordinates
          lspn=lspnsymc(isym)
          call i3minv(symlat(1,1,lspn),sym)
          s(:,:)=dble(sym(:,:))
          call r3mm(s,ainv,s)
          call r3mm(avec,s,s)
! global spin rotation of vector field
          if (ndmag.eq.3) then
! non-collinear case
            do ir=1,nrmt(is),lrstp
              do lm=1,lmmaxvr
                rvfmt(lm,ir,jas,1)=s(1,1)*rvfmt(lm,ir,ias,1) &
                                  +s(1,2)*rvfmt(lm,ir,ias,2) &
                                  +s(1,3)*rvfmt(lm,ir,ias,3)
                rvfmt(lm,ir,jas,2)=s(2,1)*rvfmt(lm,ir,ias,1) &
                                  +s(2,2)*rvfmt(lm,ir,ias,2) &
                                  +s(2,3)*rvfmt(lm,ir,ias,3)
                rvfmt(lm,ir,jas,3)=s(3,1)*rvfmt(lm,ir,ias,1) &
                                  +s(3,2)*rvfmt(lm,ir,ias,2) &
                                  +s(3,3)*rvfmt(lm,ir,ias,3)
              end do
            end do
          else
! collinear case
            do ir=1,nrmt(is),lrstp
              do lm=1,lmmaxvr
                rvfmt(lm,ir,jas,1)=s(3,3)*rvfmt(lm,ir,ias,1)
              end do
            end do
          end if
! mark atom as done
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
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
do i=1,ndmag
  call zfftifc(3,ngrid,1,zfft2(1,i))
  rvfir(:,i)=t1*dble(zfft2(:,i))
end do
deallocate(rvfmt1,rvfmt2,zfft1,zfft2)
return
end subroutine

