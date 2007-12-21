
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fxc_alda
! !INTERFACE:
subroutine kernxc
! !DESCRIPTION:
!   Computes the ALDA exchange-correlation kernel. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT. This routine is based upon the routine {\tt potxc}.
!
! !REVISION HISTORY:
!   Created March 2007 (Sagmeister)
!EOP
!BOC
  use modmain
  use modxs
  implicit none
  ! local variables
  real(8), allocatable :: dvx(:,:),dvc(:,:),rftp(:,:)
  real(8), allocatable :: f1ir(:),f2ir(:),f1mt(:,:,:),f2mt(:,:,:)
  complex(8),allocatable :: zftp(:,:)
  integer :: m,is,ia,ias,ir,itp
  allocate(rftp(lmmaxvr,1),zftp(lmmaxvr,1))
  m=max(lmmaxvr,ngrtot)
  allocate(dvx(m,1),dvc(m,1))
  !-----------------------------!
  !     muffin-tin potential    !
  !-----------------------------!
  do is=1,nspecies
     do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ir=1,nrmt(is)
           !--------------------------!
           !     spin-unpolarised     !
           !--------------------------!
           ! convert density to real space
           call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw, &
                rhomt(1,ir,ias),1,0.d0,rftp,1)
           call xcd_pwca(lmmaxvr,rftp,dvx,dvc)
           do itp=1,lmmaxvr
              rftp(itp,1)=dvx(itp,1)+dvc(itp,1)
           end do
           ! convert kernel to spherical-harmonics expansion
           zftp(:,:)=rftp(:,:)
           call zgemv('N',lmmaxvr,lmmaxvr,1.d0,zfshtvr,lmmaxvr,zftp,1,0.d0, &
                fxcmt(1,ir,ias),1)
        end do
     end do
  end do
  !--------------------------------!
  !     interstitial potential     !
  !--------------------------------!
  !--------------------------!
  !     spin-unpolarised     !
  !--------------------------!
  call xcd_pwca(ngrtot,rhoir,dvx,dvc)
  fxcir(1:ngrtot)=dvx(1:ngrtot,1)+dvc(1:ngrtot,1)
  allocate(f1ir(ngrtot),f1mt(lmmaxvr,nrmtmax,natmtot))
  allocate(f2ir(ngrtot),f2mt(lmmaxvr,nrmtmax,natmtot))
  f1ir(:)=dble(fxcir(:))
  f1mt(:,:,:)=dble(fxcmt(:,:,:))
  f2ir(:)=aimag(fxcir(:))
  f2mt(:,:,:)=aimag(fxcmt(:,:,:))
  ! symmetrise the exchange-correlation kernel
  call symrf(1,f1mt,f1ir)
  call symrf(1,f2mt,f2ir)
  ! back-substitute
  fxcir(:)=f1ir(:)+zi*f2ir(:)
  fxcmt(:,:,:)=f1mt(:,:,:)+zi*f2mt(:,:,:)
  deallocate(f1ir,f2ir,f1mt,f2mt)
  deallocate(rftp,zftp,dvx,dvc)
end subroutine kernxc
!EOC
