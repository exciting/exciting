
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genbmatk(bmt,bir,wfmt,wfir,bmat)
! calculates the magnetic field matrix elements
use modmain
implicit none
! arguments
real(8), intent(in) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
real(8), intent(in) :: bir(ngrtot,ndmag)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngrtot,nspinor,nstsv)
complex(8), intent(out) :: bmat(nstsv,nstsv)
! local variables
integer is,ia,ias,nrc,ir,irc
integer ist,jst,idm,ispn
real(8) t1
complex(8) zt1
! automatic arrays
complex(8) zflm(lmmaxvr)
! allocatable arrays
real(8), allocatable :: rvfir(:,:)
complex(8), allocatable :: zfmt(:,:,:)
complex(8), allocatable :: zfir(:,:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
allocate(rvfir(ngrtot,ndmag))
allocate(zfmt(lmmaxvr,nrcmtmax,nspinor))
allocate(zfir(ngrtot,nspinor))
! zero the matrix elements
bmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do jst=1,nstsv
  do is=1,nspecies
    nrc=nrcmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! apply magnetic field to spinor wavefunction
      do irc=1,nrc
        zfmt(:,irc,1)=bmt(:,irc,ias,ndmag)*wfmt(:,irc,ias,1,jst)
        zfmt(:,irc,2)=-bmt(:,irc,ias,ndmag)*wfmt(:,irc,ias,2,jst)
        if (ncmag) then
          zflm(:)=cmplx(bmt(:,irc,ias,1),bmt(:,irc,ias,2),8)
          zfmt(:,irc,1)=zfmt(:,irc,1)+conjg(zflm(:))*wfmt(:,irc,ias,2,jst)
          zfmt(:,irc,2)=zfmt(:,irc,2)+zflm(:)*wfmt(:,irc,ias,1,jst)
        end if
      end do
      do ist=1,jst
! compute inner product (functions are in spherical coordinates)
        do ispn=1,nspinor
          zt1=zfmtinp(.false.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
           wfmt(:,:,ias,ispn,ist),zfmt(:,:,ispn))
          bmat(ist,jst)=bmat(ist,jst)+zt1
        end do
      end do
    end do
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
do idm=1,ndmag
  rvfir(:,idm)=bir(:,idm)*cfunir(:)
end do
t1=omega/dble(ngrtot)
do jst=1,nstsv
! apply magnetic field to spinor wavefunction
  do ir=1,ngrtot
    zfir(ir,1)=rvfir(ir,ndmag)*wfir(ir,1,jst)
    zfir(ir,2)=-rvfir(ir,ndmag)*wfir(ir,2,jst)
  end do
  if (ncmag) then
    do ir=1,ngrtot
      zt1=cmplx(rvfir(ir,1),rvfir(ir,2),8)
      zfir(ir,1)=zfir(ir,1)+conjg(zt1)*wfir(ir,2,jst)
      zfir(ir,2)=zfir(ir,2)+zt1*wfir(ir,1,jst)
    end do
  end if
  do ist=1,jst
    do ispn=1,nspinor
      zt1=zdotc(ngrtot,wfir(:,ispn,ist),1,zfir(:,ispn),1)
      bmat(ist,jst)=bmat(ist,jst)+t1*zt1
    end do
  end do
end do
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    bmat(ist,jst)=conjg(bmat(jst,ist))
  end do
end do
deallocate(rvfir,zfmt,zfir)
return
end subroutine

