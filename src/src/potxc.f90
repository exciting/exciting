
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
integer is,ia,ias
integer n,nr,ir,idm,i
real(8) bext(3),t1,t2,t3,t4,ta,tb
! allocatable arrays
real(8), allocatable :: rho(:),rhoup(:),rhodn(:)
real(8), allocatable :: gvrho(:),gvup(:),gvdn(:)
real(8), allocatable :: grho(:),gup(:),gdn(:)
real(8), allocatable :: g2rho(:),g2up(:),g2dn(:)
real(8), allocatable :: g3rho(:),g3up(:),g3dn(:)
real(8), allocatable :: grho2(:),gup2(:),gdn2(:),gupdn(:)
real(8), allocatable :: ex(:),ec(:),vxc(:)
real(8), allocatable :: vx(:),vxup(:),vxdn(:)
real(8), allocatable :: vc(:),vcup(:),vcdn(:)
real(8), allocatable :: dxdg2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
real(8), allocatable :: dcdg2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
real(8), allocatable :: mag(:,:),bxc(:,:)
n=lmmaxvr*nrmtmax
allocate(rho(n),ex(n),ec(n),vxc(n))
if (associated(input%groundstate%spin)) then
  allocate(mag(n,3),bxc(n,3))
end if
n=max(n,ngrtot)
if (associated(input%groundstate%spin)) then
  allocate(rhoup(n),rhodn(n))
  allocate(vxup(n),vxdn(n),vcup(n),vcdn(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),gup(n),gdn(n))
    allocate(g2up(n),g2dn(n))
    allocate(g3rho(n),g3up(n),g3dn(n))
  else if (xcgrad.eq.2) then
    allocate(g2up(n),g2dn(n))
    allocate(gvup(3*n),gvdn(3*n))
    allocate(gup2(n),gdn2(n),gupdn(n))
    allocate(dxdgu2(n),dxdgd2(n),dxdgud(n))
    allocate(dcdgu2(n),dcdgd2(n),dcdgud(n))
  end if
else
  allocate(vx(n),vc(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),g2rho(n),g3rho(n))
  else if (xcgrad.eq.2) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
    allocate(dxdg2(n),dcdg2(n))
  end if
end if
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
do is=1,nspecies
  nr=nrmt(is)
  n=lmmaxvr*nr
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the density in spherical coordinates
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,:,ias), &
     lmmaxvr,0.d0,rho,lmmaxvr)
    if (associated(input%groundstate%spin)) then
!------------------------!
!     spin-polarised     !
!------------------------!
! magnetisation in spherical coordinates
      do idm=1,ndmag
        call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
         magmt(:,:,ias,idm),lmmaxvr,0.d0,mag(:,idm),lmmaxvr)
      end do
      if (ncmag) then
        bext(:)=input%groundstate%spin%bfieldc(:)+&
         &input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
! non-collinear (use Kubler's trick)
        do i=1,n
! compute rhoup=(rho+sgn(m.B_ext)|m|)/2 and rhodn=(rho-sgn(m.B_ext)|m|)/2
          t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
          if (xcgrad.ne.0) then
            t2=mag(i,1)*bext(1)+mag(i,2)*bext(2)+mag(i,3)*bext(3)
            if (t2.lt.0.d0) t1=-t1
          end if
          rhoup(i)=0.5d0*(rho(i)+t1)
          rhodn(i)=0.5d0*(rho(i)-t1)
        end do
      else
! collinear
        do i=1,n
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
          rhoup(i)=0.5d0*(rho(i)+mag(i,1))
          rhodn(i)=0.5d0*(rho(i)-mag(i,1))
        end do
      end if
! call the exchange-correlation interface routine
      if (xcgrad.le.0) then
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec,vxup=vxup, &
         vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.1) then
        call ggamt_sp_1(is,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
         gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
         ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.2) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
         gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
         dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2, &
         dcdgd2=dcdgd2,dcdgud=dcdgud)
        call ggamt_sp_2b(is,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
         dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
      end if
      if (ncmag) then
! non-collinear: locally spin rotate the exchange-correlation potential
        do i=1,n
          if (xctype(1).eq.100) then
            t1=vxup(i)+ec_coef*vcup(i)
            t2=vxdn(i)+ec_coef*vcdn(i)
          else
            t1=(1-ex_coef)*vxup(i)+ec_coef*vcup(i)
            t2=(1-ex_coef)*vxdn(i)+ec_coef*vcdn(i)
          end if
          vxc(i)=0.5d0*(t1+t2)
! determine the exchange-correlation magnetic field
          t3=0.5d0*(t1-t2)
          t4=rhoup(i)-rhodn(i)
          if (abs(t4).gt.1.d-8) t4=t3/t4
          bxc(i,1:3)=mag(i,1:3)*t4
        end do
      else
! collinear
        do i=1,n
          if (xctype(1).eq.100) then
            t1=vxup(i)+ec_coef*vcup(i)
            t2=vxdn(i)+ec_coef*vcdn(i)
          else
            t1=(1-ex_coef)*vxup(i)+ec_coef*vcup(i)
            t2=(1-ex_coef)*vxdn(i)+ec_coef*vcdn(i)
          end if
          vxc(i)=0.5d0*(t1+t2)
          bxc(i,1)=0.5d0*(t1-t2)
        end do
      end if
! convert field to spherical harmonics
      do idm=1,ndmag
        call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,bxc(:,idm), &
         lmmaxvr,0.d0,bxcmt(:,:,ias,idm),lmmaxvr)
      end do
      
    else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
      if (xcgrad.le.0) then
        call xcifc(xctype,n=n,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.1) then
        call ggamt_1(is,ia,grho,g2rho,g3rho)
        call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
                   ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.2) then
        call ggamt_2a(is,ia,g2rho,gvrho,grho2)
        call xcifc(xctype,n=n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
                   dxdg2=dxdg2,dcdg2=dcdg2)
        call ggamt_2b(is,g2rho,gvrho,vx,vc,dxdg2,dcdg2)
      end if
      if (xctype(1)==100) then
         vxc(1:n)=vx(1:n)+ec_coef*vc(1:n)
      else
         vxc(1:n)=(1-ex_coef)*vx(1:n)+ec_coef*vc(1:n)
      end if
    end if
    
! convert exchange and correlation energy densities to spherical harmonics
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ex,lmmaxvr, &
               0.d0,exmt(:,:,ias),lmmaxvr)
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ec,lmmaxvr, &
               0.d0,ecmt(:,:,ias),lmmaxvr)
               
    if (xctype(1).ne.100) exmt(:,:,ias) = (1.d0-ex_coef)*exmt(:,:,ias)
    ecmt(:,:,ias) = ec_coef*ecmt(:,:,ias)
    
! convert exchange-correlation potential to spherical harmonics
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,vxc,lmmaxvr, &
                0.d0,vxcmt(:,:,ias),lmmaxvr)
      
  end do
end do

!------------------------------------------!
!     interstitial potential and field     !
!------------------------------------------!

if (associated(input%groundstate%spin)) then
  !------------------------!
  !     spin-polarised     !
  !------------------------!
  if (ncmag) then
    ! non-collinear
    do ir=1,ngrtot
      t1=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
      if (xcgrad.ne.0) then
        t2 = magir(ir,1)*input%groundstate%spin%bfieldc(1)+ &
        &    magir(ir,2)*input%groundstate%spin%bfieldc(2)+ &
        &    magir(ir,3)*input%groundstate%spin%bfieldc(3)
        if (t2.lt.0.d0) t1=-t1
      end if
      rhoup(ir)=0.5d0*(rhoir(ir)+t1)
      rhodn(ir)=0.5d0*(rhoir(ir)-t1)
    end do
  else
    ! collinear
    do ir=1,ngrtot
      rhoup(ir)=0.5d0*(rhoir(ir)+magir(ir,1))
      rhodn(ir)=0.5d0*(rhoir(ir)-magir(ir,1))
    end do
  end if
  
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,ex=exir,ec=ecir, &
    &          vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.1) then
    call ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
    &          gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
    &          ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
  else if (xcgrad.eq.2) then
    call ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
    call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
    &          gupdn=gupdn,ex=exir,ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
    &          dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
    &          dcdgud=dcdgud)
    call ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
    &                dxdgud,dcdgu2,dcdgd2,dcdgud)
  end if
  
  if (ncmag) then
    ! non-collinear: spin rotate the local exchange potential
    do ir = 1, ngrtot
      if (xctype(1).eq.100) then
        t1 = vxup(ir)+ec_coef*vcup(ir)
        t2 = vxdn(ir)+ec_coef*vcdn(ir)
      else
        t1 = (1.d0-ex_coef)*vxup(ir)+ec_coef*vcup(ir)
        t2 = (1.d0-ex_coef)*vxdn(ir)+ec_coef*vcdn(ir)
      end if
      vxcir(ir) = 0.5d0*(t1+t2)
      ! determine the exchange-correlation magnetic field
      t3 = 0.5d0*(t1-t2)
      t4 = rhoup(ir)-rhodn(ir)
      if (abs(t4).gt.1.d-8) t4 = t3/t4
      bxcir(ir,:) = magir(ir,:)*t4
    end do
  else
    ! collinear
    do ir = 1, ngrtot
      if (xctype(1).eq.100) then
        t1 = vxup(ir)+ec_coef*vcup(ir)
        t2 = vxdn(ir)+ec_coef*vcdn(ir)
      else
        t1 = (1.d0-ex_coef)*vxup(ir)+ec_coef*vcup(ir)
        t2 = (1.d0-ex_coef)*vxdn(ir)+ec_coef*vcdn(ir)
      end if
      vxcir(ir) = 0.5d0*(t1+t2)
      bxcir(ir,1) = 0.5d0*(t1-t2)
    end do
  end if
  if (xctype(1).ne.100) exir(:) = (1.d0-ex_coef)*exir(:)
  ecir(:) = ec_coef*ecir(:)
  
else

  !--------------------------!
  !     spin-unpolarised     !
  !--------------------------!
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngrtot,rho=rhoir,ex=exir,ec=ecir,vx=vx,vc=vc)
  else if (xcgrad.eq.1) then
    call ggair_1(grho,g2rho,g3rho)
  call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, &
  &          ex=exir,ec=ecir,vx=vx,vc=vc)
  else if (xcgrad.eq.2) then
    call ggair_2a(g2rho,gvrho,grho2)
    call xcifc(xctype,n=ngrtot,rho=rhoir,grho2=grho2,ex=exir,ec=ecir,vx=vx, &
    &          vc=vc,dxdg2=dxdg2,dcdg2=dcdg2)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdg2,dcdg2)
  end if
  if (xctype(1).ne.100) then
    vxcir(1:ngrtot) = (1.d0-ex_coef)*vx(1:ngrtot)+ec_coef*vc(1:ngrtot)
    exir(:) = (1.d0-ex_coef)*exir(:)
  else
    vxcir(1:ngrtot) = vx(1:ngrtot)+ec_coef*vc(1:ngrtot)
  end if
  ecir(:) = ec_coef*ecir(:)

end if ! spin case

! OEP - EXX/HYBRIDS
If  (associated(input%groundstate%OEP)) call oepmain

!-----------------------------------------------
! symmetrise the exchange-correlation potential
!-----------------------------------------------
call symrf(1,vxcmt,vxcir)
if (associated(input%groundstate%spin)) then
  ! remove the source contribution if required
  if (input%groundstate%nosource) call projsbf
  ! symmetrise the exchange-correlation effective field
  call symrvf(1,bxcmt,bxcir)
end if

! Clear memory
deallocate(rho,ex,ec,vxc)
if (associated(input%groundstate%spin)) then
  deallocate(mag,bxc)
  deallocate(rhoup,rhodn,vxup,vxdn,vcup,vcdn)
  if (xcgrad.eq.1) then
    deallocate(grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
  else if (xcgrad.eq.2) then
    deallocate(g2up,g2dn)
    deallocate(gvup,gvdn)
    deallocate(gup2,gdn2,gupdn)
    deallocate(dxdgu2,dxdgd2,dxdgud)
    deallocate(dcdgu2,dcdgd2,dcdgud)
  end if
else
  deallocate(vx,vc)
  if (xcgrad.eq.1) then
    deallocate(grho,g2rho,g3rho)
  else if (xcgrad.eq.2) then
    deallocate(g2rho,gvrho,grho2)
    deallocate(dxdg2,dcdg2)
  end if
end if

return
end subroutine
!EOC

