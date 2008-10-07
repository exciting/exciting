
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potxc
! !INTERFACE:
subroutine potxc
! !USES:
use modmain
use modxcifc
! !DESCRIPTION:
!   Computes the exchange-correlation potential and energy density. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer m,is,ia,ias,ir,itp,i
real(8) t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: rftp1(:,:)
real(8), allocatable :: rftp2(:,:)
real(8), allocatable :: rfir(:,:)
real(8), allocatable :: ex(:)
real(8), allocatable :: ec(:)
real(8), allocatable :: vx(:,:)
real(8), allocatable :: vc(:,:)
real(8), allocatable :: grhomt(:,:)
real(8), allocatable :: gupmt(:,:)
real(8), allocatable :: gdnmt(:,:)
real(8), allocatable :: g2upmt(:,:)
real(8), allocatable :: g2dnmt(:,:)
real(8), allocatable :: g3rhomt(:,:)
real(8), allocatable :: g3upmt(:,:)
real(8), allocatable :: g3dnmt(:,:)
real(8), allocatable :: grhoir(:)
real(8), allocatable :: gupir(:)
real(8), allocatable :: gdnir(:)
real(8), allocatable :: g2upir(:)
real(8), allocatable :: g2dnir(:)
real(8), allocatable :: g3rhoir(:)
real(8), allocatable :: g3upir(:)
real(8), allocatable :: g3dnir(:)
allocate(ex(lmmaxvr))
allocate(ec(lmmaxvr))
allocate(rftp1(lmmaxvr,4))
allocate(rftp2(lmmaxvr,2))
m=max(lmmaxvr,ngrtot)
if (spinpol) then
  allocate(rfir(ngrtot,2))
  allocate(vx(m,2),vc(m,2))
else
  allocate(vx(m,1),vc(m,1))
end if
if (xcgrad.eq.1) then
  allocate(grhomt(lmmaxvr,nrmtmax))
  allocate(gupmt(lmmaxvr,nrmtmax))
  allocate(gdnmt(lmmaxvr,nrmtmax))
  allocate(g2upmt(lmmaxvr,nrmtmax))
  allocate(g2dnmt(lmmaxvr,nrmtmax))
  allocate(g3rhomt(lmmaxvr,nrmtmax))
  allocate(g3upmt(lmmaxvr,nrmtmax))
  allocate(g3dnmt(lmmaxvr,nrmtmax))
  allocate(grhoir(ngrtot))
  allocate(gupir(ngrtot))
  allocate(gdnir(ngrtot))
  allocate(g2upir(ngrtot))
  allocate(g2dnir(ngrtot))
  allocate(g3rhoir(ngrtot))
  allocate(g3upir(ngrtot))
  allocate(g3dnir(ngrtot))
end if
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! gradients for GGA if required
    if (xcgrad.eq.1) call ggamt(is,ia,grhomt,gupmt,gdnmt,g2upmt,g2dnmt, &
     g3rhomt,g3upmt,g3dnmt)
    do ir=1,nrmt(is)
      if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
        if (ncmag) then
! non-collinear
          do i=1,3
            call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
             magmt(:,ir,ias,i),1,0.d0,rftp1(:,i),1)
          end do
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,ir,ias), &
           1,0.d0,rftp1(:,4),1)
! compute (rho+|m|)/2 and (rho-|m|)/2
          do itp=1,lmmaxvr
            t1=rftp1(itp,4)
            t2=sqrt(rftp1(itp,1)**2+rftp1(itp,2)**2+rftp1(itp,3)**2)
            rftp2(itp,1)=0.5d0*(t1+t2)
            rftp2(itp,2)=0.5d0*(t1-t2)
          end do
        else
! collinear
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,ir,ias), &
           1,0.d0,rftp1(:,1),1)
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
           magmt(:,ir,ias,1),1,0.d0,rftp1(:,2),1)
          do itp=1,lmmaxvr
            t1=rftp1(itp,1)
            t2=rftp1(itp,2)
            rftp2(itp,1)=0.5d0*(t1+t2)
            rftp2(itp,2)=0.5d0*(t1-t2)
          end do
        end if
        if (xcgrad.le.0) then
          call xcifc(xctype,n=lmmaxvr,rhoup=rftp2(:,1),rhodn=rftp2(:,2),ex=ex, &
           ec=ec,vxup=vx(:,1),vxdn=vx(:,2),vcup=vc(:,1),vcdn=vc(:,2))
        else
          call xcifc(xctype,n=lmmaxvr,rhoup=rftp2(:,1),rhodn=rftp2(:,2),ex=ex, &
           ec=ec,vxup=vx(:,1),vxdn=vx(:,2),vcup=vc(:,1),vcdn=vc(:,2), &
           grho=grhomt(:,ir),gup=gupmt(:,ir),gdn=gdnmt(:,ir), &
           g2up=g2upmt(:,ir),g2dn=g2dnmt(:,ir),g3rho=g3rhomt(:,ir), &
           g3up=g3upmt(:,ir),g3dn=g3dnmt(:,ir))
        end if
        if (ncmag) then
! non-collinear: spin rotate the local exchange-correlation potential
          do itp=1,lmmaxvr
            t1=vx(itp,1)+vc(itp,1)
            t2=vx(itp,2)+vc(itp,2)
            t3=0.5d0*(t1-t2)
! determine |m| again
            t4=2.d0*rftp2(itp,1)-rftp1(itp,4)
            if (t4.gt.1.d-8) t4=t3/t4
            rftp1(itp,1:3)=rftp1(itp,1:3)*t4
            rftp1(itp,4)=0.5d0*(t1+t2)
          end do
          do i=1,3
            call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp1(:,i),1, &
             0.d0,bxcmt(:,ir,ias,i),1)
          end do
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp1(:,4),1, &
           0.d0,vxcmt(:,ir,ias),1)
        else
! collinear
          do itp=1,lmmaxvr
            t1=vx(itp,1)+vc(itp,1)
            t2=vx(itp,2)+vc(itp,2)
            rftp1(itp,1)=0.5d0*(t1+t2)
            rftp1(itp,2)=0.5d0*(t1-t2)
          end do
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp1(:,1),1, &
           0.d0,vxcmt(:,ir,ias),1)
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp1(:,2),1, &
           0.d0,bxcmt(:,ir,ias,1),1)
        end if
      else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,ir,ias), &
         1,0.d0,rftp1,1)
        if (xcgrad.le.0) then
          call xcifc(xctype,n=lmmaxvr,rho=rftp1,ex=ex,ec=ec,vx=vx,vc=vc)
        else
          call xcifc(xctype,n=lmmaxvr,rho=rftp1,ex=ex,ec=ec,vx=vx,vc=vc, &
           grho=gupmt(:,ir),g2rho=g2upmt(:,ir),g3rho=g3upmt(:,ir))
        end if
        do itp=1,lmmaxvr
          rftp1(itp,1)=vx(itp,1)+vc(itp,1)
        end do
! exchange-correlation potential
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp1,1,0.d0, &
         vxcmt(:,ir,ias),1)
      end if
! convert energy densities from spherical coordinates to spherical harmonics
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ex,1,0.d0, &
       exmt(:,ir,ias),1)
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ec,1,0.d0, &
       ecmt(:,ir,ias),1)
    end do
  end do
end do
!--------------------------------!
!     interstitial potential     !
!--------------------------------!
! gradients for GGA if required
if (xcgrad.eq.1) call ggair(grhoir,gupir,gdnir,g2upir,g2dnir,g3rhoir,g3upir, &
 g3dnir)
if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
  do ir=1,ngrtot
    if (ncmag) then
! non-collinear
      t1=rhoir(ir)
! compute |m|
      t2=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
    else
! collinear
      t1=rhoir(ir)
      t2=magir(ir,1)
    end if
    rfir(ir,1)=0.5d0*(t1+t2)
    rfir(ir,2)=0.5d0*(t1-t2)
  end do
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngrtot,rhoup=rfir(:,1),rhodn=rfir(:,2),ex=exir, &
     ec=ecir,vxup=vx(:,1),vxdn=vx(:,2),vcup=vc(:,1),vcdn=vc(:,2))
  else
    call xcifc(xctype,n=ngrtot,rhoup=rfir(:,1),rhodn=rfir(:,2),ex=exir, &
     ec=ecir,vxup=vx(:,1),vxdn=vx(:,2),vcup=vc(:,1),vcdn=vc(:,2),grho=grhoir, &
     gup=gupir,gdn=gdnir,g2up=g2upir,g2dn=g2dnir,g3rho=g3rhoir,g3up=g3upir, &
     g3dn=g3dnir)
  end if
  if (ncmag) then
! non-collinear: spin rotate the local exchange potential
    do ir=1,ngrtot
      t1=vx(ir,1)+vc(ir,1)
      t2=vx(ir,2)+vc(ir,2)
      t3=0.5d0*(t1-t2)
! determine |m| again
      t4=2.d0*rfir(ir,1)-rhoir(ir)
      if (t4.gt.1.d-8) t4=t3/t4
      bxcir(ir,:)=magir(ir,:)*t4
      vxcir(ir)=0.5d0*(t1+t2)
    end do
  else
! collinear
    do ir=1,ngrtot
      t1=vx(ir,1)+vc(ir,1)
      t2=vx(ir,2)+vc(ir,2)
      vxcir(ir)=0.5d0*(t1+t2)
      bxcir(ir,1)=0.5d0*(t1-t2)
    end do
  end if
else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngrtot,rho=rhoir,ex=exir,ec=ecir,vx=vx,vc=vc)
  else
    call xcifc(xctype,n=ngrtot,rho=rhoir,ex=exir,ec=ecir,vx=vx,vc=vc, &
     grho=gupir,g2rho=g2upir,g3rho=g3upir)
  end if
  vxcir(1:ngrtot)=vx(1:ngrtot,1)+vc(1:ngrtot,1)
end if
! optimised effective potential
if (xctype.lt.0) call oepmain
! symmetrise the exchange-correlation potential
call symrf(1,vxcmt,vxcir)
if (spinpol) then
! remove the source contribution if required
  if (nosource) call projsbf
! symmetrise the exchange-correlation effective field
  call symrvf(1,bxcmt,bxcir)
end if
deallocate(ex,ec,vx,vc,rftp1,rftp2)
if (spinpol) deallocate(rfir)
if (xcgrad.eq.1) then
  deallocate(grhomt,gupmt,gdnmt,g2upmt,g2dnmt,g3rhomt,g3upmt,g3dnmt)
  deallocate(grhoir,gupir,gdnir,g2upir,g2dnir,g3rhoir,g3upir,g3dnir)
end if
return
end subroutine
!EOC

