
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fxc_alda_check
  !
  ! For checking the validity of the analytical expression for the
  ! ALDA exchange-correlation kernel
  ! Apr. 2007 (Sagmeister)
  !
  use modmain
  use modxcifc
  use modxs
  use modfxcifc
  implicit none
  real(8), allocatable :: vx(:),ex(:)
  real(8), allocatable :: vc(:),ec(:)
  real(8), allocatable :: dvx(:),dvc(:),vxc(:),dvxc2(:),dvxc(:),cf(:,:)
  integer :: m,igr,nrho,irh
  real(8), allocatable :: rhogr(:)
  real(8) :: rhoint(2)

  call init0
  call init1
  call readstate
  call init2xs

  m=max(lmmaxvr,ngrtot)

  nrho=100000
!!$  rhoint(1)=1.0d-3
!!$  rhoint(2)=1.d0
  rhoint(1)=1.d0
  rhoint(2)=7400.0d0
  
  
  allocate(rhogr(nrho))
  allocate(dvx(nrho),dvc(nrho),dvxc2(nrho),dvxc(nrho),cf(3,nrho))
  allocate(ex(nrho),ec(nrho),vx(nrho),vc(nrho),vxc(nrho))

  do irh=1,nrho
     rhogr(irh)=rhoint(1)+(rhoint(2)-rhoint(1))*dble(irh-1)/dble(nrho)
  end do

  call xcifc(xctype,n=nrho,rho=rhogr,ex=ex,ec=ec,vx=vx,vc=vc)
  vxc(:)=vx(:)+vc(:)

  ! use analytic expression of xc-kernel
  call xcd_pwca(nrho,rhogr,dvx,dvc)
  dvxc(:)=dvx(:)+dvc(:)

  ! numerical differentiation
  call fderiv(1,nrho,rhogr,vxc,dvxc2,cf)

  ! plot both versions

  do irh=1,nrho
     write(2000,'(i6,100g18.10)') irh,rhogr(irh),vxc(irh),dvxc2(irh),dvxc(irh)
  end do

  write(*,*) 'minimum interstital density:',minval(rhoir)
  write(*,*) 'minimum muffin-tin density :',minval(rhomt)
  write(*,*) 'maximum interstital density:',maxval(rhoir)
  write(*,*) 'maximum muffin-tin density :',maxval(rhomt)
  write(*,*)
  write(*,*) 'minimum interstital potential:',minval(vxcir)
  write(*,*) 'minimum muffin-tin potential :',minval(vxcmt)
  write(*,*) 'maximum interstital potential:',maxval(vxcir)
  write(*,*) 'maximum muffin-tin potential :',maxval(vxcmt)


end subroutine fxc_alda_check
