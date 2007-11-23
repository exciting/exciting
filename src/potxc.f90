
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
real(8) w1,w2,t1,t2,t3
complex(8) a(2,2)
! allocatable arrays
real(8), allocatable :: rflm(:,:)
real(8), allocatable :: rftp(:,:)
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
real(8), allocatable :: cs1(:)
complex(8), allocatable :: sn1(:)
allocate(ex(lmmaxvr))
allocate(ec(lmmaxvr))
m=max(lmmaxvr,ngrtot)
if (spinpol) then
  allocate(rflm(lmmaxvr,4))
  allocate(rftp(lmmaxvr,4))
  allocate(rfir(ngrtot,2))
  allocate(vx(m,2),vc(m,2))
  allocate(cs1(m),sn1(m))
else
  allocate(rftp(lmmaxvr,1))
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
        if (ndmag.eq.3) then
! non-collinear
          rflm(:,1)=0.5d0*(rhomt(:,ir,ias)+magmt(:,ir,ias,3))
          rflm(:,2)=0.5d0*(rhomt(:,ir,ias)-magmt(:,ir,ias,3))
          rflm(:,3)=0.5d0*magmt(:,ir,ias,1)
          rflm(:,4)=-0.5d0*magmt(:,ir,ias,2)
          do i=1,4
            call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,rflm(1,i),1, &
             0.d0,rftp(1,i),1)
          end do
! diagonalise the local spin density
          do itp=1,lmmaxvr
            a(1,1)=rftp(itp,1)
            a(2,2)=rftp(itp,2)
            a(1,2)=cmplx(rftp(itp,3),rftp(itp,4),8)
            call zlaev2(a(1,1),a(1,2),a(2,2),rftp(itp,1),rftp(itp,2),cs1(itp), &
             sn1(itp))
          end do
        else
! collinear
          rflm(:,1)=0.5d0*(rhomt(:,ir,ias)+magmt(:,ir,ias,1))
          rflm(:,2)=0.5d0*(rhomt(:,ir,ias)-magmt(:,ir,ias,1))
          do i=1,2
            call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,rflm(1,i),1, &
             0.d0,rftp(1,i),1)
          end do
        end if
        if (xcgrad.le.0) then
          call xcifc(xctype,n=lmmaxvr,rhoup=rftp(1,1),rhodn=rftp(1,2),ex=ex, &
           ec=ec,vxup=vx(1,1),vxdn=vx(1,2),vcup=vc(1,1),vcdn=vc(1,2))
        else
          call xcifc(xctype,n=lmmaxvr,rhoup=rftp(1,1),rhodn=rftp(1,2),ex=ex, &
           ec=ec,vxup=vx(1,1),vxdn=vx(1,2),vcup=vc(1,1),vcdn=vc(1,2), &
           grho=grhomt(1,ir),gup=gupmt(1,ir),gdn=gdnmt(1,ir), &
           g2up=g2upmt(1,ir),g2dn=g2dnmt(1,ir),g3rho=g3rhomt(1,ir), &
           g3up=g3upmt(1,ir),g3dn=g3dnmt(1,ir))
        end if
        if (ndmag.eq.3) then
! non-collinear: spin rotate the local exchange-correlation potential
          do itp=1,lmmaxvr
            w1=vx(itp,1)+vc(itp,1)
            w2=vx(itp,2)+vc(itp,2)
            t1=cs1(itp)
            t2=dble(sn1(itp))
            t3=aimag(sn1(itp))
            rftp(itp,1)=t1*t2*(w1-w2)
            rftp(itp,2)=t1*t3*(w1-w2)
            rftp(itp,3)=0.5d0*(t1**2-(t2**2+t3**2))*(w1-w2)
            rftp(itp,4)=0.5d0*(w1+w2)
          end do
          do i=1,3
            call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp(1,i),1, &
             0.d0,bxcmt(1,ir,ias,i),1)
          end do
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp(1,4),1, &
           0.d0,vxcmt(1,ir,ias),1)
        else
! collinear
          do itp=1,lmmaxvr
            rftp(itp,1)=0.5d0*(vx(itp,1)+vc(itp,1)+vx(itp,2)+vc(itp,2))
            rftp(itp,2)=0.5d0*(vx(itp,1)+vc(itp,1)-vx(itp,2)-vc(itp,2))
          end do
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp(1,1),1, &
           0.d0,vxcmt(1,ir,ias),1)
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp(1,2),1, &
           0.d0,bxcmt(1,ir,ias,1),1)
        end if
      else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,rhomt(1,ir,ias), &
         1,0.d0,rftp,1)
        if (xcgrad.le.0) then
          call xcifc(xctype,n=lmmaxvr,rho=rftp,ex=ex,ec=ec,vx=vx,vc=vc)
        else
          call xcifc(xctype,n=lmmaxvr,rho=rftp,ex=ex,ec=ec,vx=vx,vc=vc, &
           grho=gupmt(1,ir),g2rho=g2upmt(1,ir),g3rho=g3upmt(1,ir))
        end if
        do itp=1,lmmaxvr
          rftp(itp,1)=vx(itp,1)+vc(itp,1)
        end do
! exchange-correlation potential
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp,1,0.d0, &
         vxcmt(1,ir,ias),1)
      end if
! convert energy densities from spherical coordinates to spherical harmonics
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ex,1,0.d0, &
       exmt(1,ir,ias),1)
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ec,1,0.d0, &
       ecmt(1,ir,ias),1)
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
    if (ndmag.eq.3) then
! non-collinear
      a(1,1)=0.5d0*(rhoir(ir)+magir(ir,3))
      a(2,2)=0.5d0*(rhoir(ir)-magir(ir,3))
      a(1,2)=cmplx(0.5d0*magir(ir,1),-0.5d0*magir(ir,2),8)
! diagonalise the local spin density
      call zlaev2(a(1,1),a(1,2),a(2,2),rfir(ir,1),rfir(ir,2),cs1(ir),sn1(ir))
    else
! collinear
      rfir(ir,1)=0.5d0*(rhoir(ir)+magir(ir,1))
      rfir(ir,2)=0.5d0*(rhoir(ir)-magir(ir,1))
    end if
  end do
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngrtot,rhoup=rfir(1,1),rhodn=rfir(1,2),ex=exir, &
     ec=ecir,vxup=vx(1,1),vxdn=vx(1,2),vcup=vc(1,1),vcdn=vc(1,2))
  else
    call xcifc(xctype,n=ngrtot,rhoup=rfir(1,1),rhodn=rfir(1,2),ex=exir, &
     ec=ecir,vxup=vx(1,1),vxdn=vx(1,2),vcup=vc(1,1),vcdn=vc(1,2),grho=grhoir, &
     gup=gupir,gdn=gdnir,g2up=g2upir,g2dn=g2dnir,g3rho=g3rhoir,g3up=g3upir, &
     g3dn=g3dnir)
  end if
  if (ndmag.eq.3) then
! non-collinear: spin rotate the local exchange potential
    do ir=1,ngrtot
      w1=vx(ir,1)+vc(ir,1)
      w2=vx(ir,2)+vc(ir,2)
      t1=cs1(ir)
      t2=dble(sn1(ir))
      t3=aimag(sn1(ir))
      bxcir(ir,1)=t1*t2*(w1-w2)
      bxcir(ir,2)=t1*t3*(w1-w2)
      bxcir(ir,3)=0.5d0*(t1**2-(t2**2+t3**2))*(w1-w2)
      vxcir(ir)=0.5d0*(w1+w2)
    end do
  else
! collinear
    do ir=1,ngrtot
      vxcir(ir)=0.5d0*(vx(ir,1)+vc(ir,1)+vx(ir,2)+vc(ir,2))
      bxcir(ir,1)=0.5d0*(vx(ir,1)+vc(ir,1)-vx(ir,2)-vc(ir,2))
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
deallocate(ex,ec,vx,vc,rftp)
if (spinpol) deallocate(rflm,rfir,cs1,sn1)
if (xcgrad.eq.1) then
  deallocate(grhomt,gupmt,gdnmt,g2upmt,g2dnmt,g3rhomt,g3upmt,g3dnmt)
  deallocate(grhoir,gupir,gdnir,g2upir,g2dnir,g3rhoir,g3upir,g3dnir)
end if
return
end subroutine
!EOC

