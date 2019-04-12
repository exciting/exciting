
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
integer n,nr,ir,idm,i,j
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
real(8), allocatable :: exsr(:),vxsr(:),vxsrup(:),vxsrdn(:),v2xsr(:),v2xsrup(:),v2xsrdn(:)!,gv2xsr(:)
real(8), allocatable :: vc(:),vcup(:),vcdn(:)
!short-range energy and potential needed for hybrids (HSE)
!real(8), allocatable :: vxsr(:), exsr(:)
real(8), allocatable :: dxdg2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
real(8), allocatable :: dcdg2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
real(8), allocatable :: mag(:,:),bxc(:,:)

n=lmmaxvr*nrmtmax
!!CECI allocate in the following part also exsr(n) and vxsr(n), and add in the rest the case for HSE sr part
!if (allocated(rho)) deallocate(rho)
!if (allocated(ex)) deallocate(ex)
!if (allocated(ec)) deallocate(ec)
!if (allocated(vxc)) deallocate(vxc)
allocate(rho(n),ex(n),ec(n),vxc(n))
if (associated(input%groundstate%spin)) then
  allocate(mag(n,3),bxc(n,3))
end if
n=max(n,ngrtot)
if (associated(input%groundstate%spin)) then
  allocate(rhoup(n),rhodn(n))
  allocate(vxup(n),vxdn(n),vcup(n),vcdn(n))
  !allocate(vxsrup(n),vxsrdn(n),v2xsrup(n),v2xsrdn(n))
  if (xcgrad.eq.1) then
    allocate(grho(n),gup(n),gdn(n))
    allocate(g2up(n),g2dn(n))
    allocate(g3rho(n),g3up(n),g3dn(n))
    !!CECI for sure not correct but we need to add
    if (xctype(1)==23.or.xctype(1)==408) then
       if (allocated(exsr)) deallocate(exsr)
       if (allocated(vxsr)) deallocate(vxsr)
       if (allocated(vxsrup)) deallocate(vxsrup)
       if (allocated(vxsrdn)) deallocate(vxsrdn)
       if (allocated(v2xsr)) deallocate(v2xsr)
       if (allocated(v2xsrup)) deallocate(v2xsrup)
       if (allocated(v2xsrdn)) deallocate(v2xsrdn)
       allocate(exsr(n),vxsr(n),vxsrup(n),vxsrdn(n),v2xsr(n),v2xsrup(n),v2xsrdn(n))!,gv2xsr(n))
    endif
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
    !if (xctype(1)==408) then
   !    allocate(vxsr(n),exsr(n))
   ! endif
    if (xctype(1)==23.or.xctype(1)==408) then
       if (allocated(exsr)) deallocate(exsr)
       if (allocated(vxsr)) deallocate(vxsr)
       if (allocated(v2xsr)) deallocate(v2xsr)
       allocate(exsr(n),vxsr(n),v2xsr(n))!,gv2xsr(n))
    endif
  else if (xcgrad.eq.2) then
    allocate(g2rho(n),gvrho(3*n),grho2(n))
    allocate(dxdg2(n),dcdg2(n))
  end if
end if

 
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
write(*,*) "Hello21"
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
write(*,*) "Hello22"
! magnetisation in spherical coordinates
      do idm=1,ndmag
        call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
         magmt(:,:,ias,idm),lmmaxvr,0.d0,mag(:,idm),lmmaxvr)
      end do
      if (ncmag) then
        bext(:)=input%groundstate%spin%bfieldc(:)+&
         &input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
write(*,*) "Hello23"
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
write(*,*) "Hello24"
! collinear
        do i=1,n
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
          rhoup(i)=0.5d0*(rho(i)+mag(i,1))
          rhodn(i)=0.5d0*(rho(i)-mag(i,1))
        end do
      end if
write(*,*) "Hello25"
! call the exchange-correlation interface routine
      if (xcgrad.le.0) then
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec,vxup=vxup, &
         vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.1) then
      !CECI PROBABLY YOU HAVE TO ADD HERE
        call ggamt_sp_1(is,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
        if (xctype(1)==23) then
write(*,*) "Hello25"
           call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
             gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
             ec=ec,exsr=exsr,vxsrup=vxsrup,vxsrdn=vxsrdn,v2xsrup=v2xsrup, v2xsrdn=v2xsrdn,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
           !call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
           !        ec=ec,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
write(*,*) "Hello25"
           call gv2xmt_spin(is,ia,gup,vxsrup,v2xsrup)
           !call gv2xmt(is,ia,grho,vxsrup,v2xsrup)
write(*,*) "Hello25"
           call gv2xmt_spin(is,ia,gdn,vxsrdn,v2xsrdn)
           !call gv2xmt(is,ia,grho,vxsrdn,v2xsrdn)
           vxup=vxsrup
           vxdn=vxsrdn
           ex=exsr
        elseif (xctype(1)==408) then
write(*,*) "Hello26"
           call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
             gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
             ec=ec,exsr=exsr,vxsrup=vxsrup,vxsrdn=vxsrdn,v2xsrup=v2xsrup, v2xsrdn=v2xsrdn,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
write(*,*) "Hello27"
           call gv2xmt_spin(is,ia,gup,vxsrup,v2xsrup)
           call gv2xmt_spin(is,ia,gdn,vxsrdn,v2xsrdn)
write(*,*) "Hello28"
           !call gv2xmt_sp(is,ia,gup,,gdn,vxsrup,vxsrdn,v2xsrup,v2xsrdn)
           !call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
           !        ec=ec,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
           !call gv2xmt(is,ia,grho,vxsr,v2xsr)
        else
           call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
            gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
            ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
        !  call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
        !           ec=ec,vx=vx,vc=vc)
        endif
        !call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
        ! gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex, &
        ! ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      else if (xcgrad.eq.2) then
        call ggamt_sp_2a(is,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
        call xcifc(xctype,n=n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
         gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
         dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2, &
         dcdgd2=dcdgd2,dcdgud=dcdgud)
        call ggamt_sp_2b(is,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
         dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
      end if
write(*,*) "Hello26"
      if (ncmag) then
! non-collinear: locally spin rotate the exchange-correlation potential
        !!CECI ALSO HERE
        do i=1,n
          if (xctype(1).eq.100) then
            t1=vxup(i)+ec_coef*vcup(i)
            t2=vxdn(i)+ec_coef*vcdn(i)
          elseif (xctype(1)==408) then !HSE
            t1=vcup(i)+vxup(i)-ex_coef*vxsrup(i)
            t2=vcdn(i)+vxdn(i)-ex_coef*vxsrdn(i)
            !vxc(1:n)=vc(1:n)+vx(1:n)-ex_coef*vxsr(1:n)     
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
write(*,*) "Hello27"
! collinear
        do i=1,n
          if (xctype(1).eq.100) then
            t1=vxup(i)+ec_coef*vcup(i)
            t2=vxdn(i)+ec_coef*vcdn(i)
          elseif (xctype(1)==408) then !HSE
            t1=vcup(i)+vxup(i)-ex_coef*vxsrup(i)
            t2=vcdn(i)+vxdn(i)-ex_coef*vxsrdn(i)
          else
            t1=(1-ex_coef)*vxup(i)+ec_coef*vcup(i)
            t2=(1-ex_coef)*vxdn(i)+ec_coef*vcdn(i)
          end if
          vxc(i)=0.5d0*(t1+t2)
          bxc(i,1)=0.5d0*(t1-t2)
        end do
      end if
write(*,*) "Hello27"
! convert field to spherical harmonics
      do idm=1,ndmag
        call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,bxc(:,idm), &
         lmmaxvr,0.d0,bxcmt(:,:,ias,idm),lmmaxvr)
      end do
write(*,*) "Hello28"
    else
!--------------------------!
!     spin-unpolarised     !
!-------------------------!
write(*,*) "Hello29"
      if (xcgrad.le.0) then
        call xcifc(xctype,n=n,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
      else if (xcgrad.eq.1) then
        call ggamt_1(is,ia,grho,g2rho,g3rho)
        !if (xctype(1)==408) then
   !       call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho, ex=ex,&
   !           ec=ec,exsr=exsr, vx=vx, vc=vc, vxsr=vxsr)
   !     else 
        if (xctype(1)==23) then
           call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
                   ec=ec,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
           call gv2xmt(is,ia,grho,vxsr,v2xsr)
           vx=vxsr
           ex=exsr 
        elseif (xctype(1)==408) then
           call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
                   ec=ec,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
           call gv2xmt(is,ia,grho,vxsr,v2xsr)
        else
           call xcifc(xctype,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
                   ec=ec,vx=vx,vc=vc)
        endif
      else if (xcgrad.eq.2) then
        call ggamt_2a(is,ia,g2rho,gvrho,grho2)
        call xcifc(xctype,n=n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
                   dxdg2=dxdg2,dcdg2=dcdg2)
        call ggamt_2b(is,g2rho,gvrho,vx,vc,dxdg2,dcdg2)
      end if
      if (xctype(1)==100) then
         vxc(1:n)=vx(1:n)+ec_coef*vc(1:n)
      elseif (xctype(1)==408) then !HSE
         vxc(1:n)=vc(1:n)+vx(1:n)-ex_coef*vxsr(1:n)     
      else
         vxc(1:n)=(1-ex_coef)*vx(1:n)+ec_coef*vc(1:n)
      end if
    end if


!CECIIII ebfi of the first part
write(*,*) "Hello29"    
    if ((xctype(1).ne.100).and.(xctype(1).ne.408)) then 
       ex(1:n) = (1.d0-ex_coef)*ex(1:n)
write(*,*) "Hello29"
    elseif (xctype(1)==408) then !HSE
       ex(1:n) = ex(1:n)-ex_coef*exsr(1:n)
    end if
    ec(1:n) = ec_coef*ec(1:n)
! convert exchange and correlation energy densities to spherical harmonics
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ex,lmmaxvr, &
               0.d0,exmt(:,:,ias),lmmaxvr)
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,ec,lmmaxvr, &
               0.d0,ecmt(:,:,ias),lmmaxvr)
!! CECI LOOK HERE               
!  if ((xctype(1).ne.100)) then !CECItest
!  if ((xctype(1).ne.100).and.(xctype(1).ne.408)) then
!     exmt(:,:,ias) = (1.d0-ex_coef)*exmt(:,:,ias)
!    !CECI:test
!  elseif (xctype(1)==408) then !HSE
!     call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,exsr,lmmaxvr, &
!               0.d0,exmtsr(:,:,ias),lmmaxvr)
!     exmt(:,:,ias) = exmt(:,:,ias)-ex_coef*exsrmt(:,:,ias)
!  end if
!    ecmt(:,:,ias) = ec_coef*ecmt(:,:,ias)
    
! convert exchange-correlation potential to spherical harmonics
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,vxc,lmmaxvr, &
                0.d0,vxcmt(:,:,ias),lmmaxvr)
  
      
  end do
end do
!------------------------------------------!
!     interstitial potential and field     !
!------------------------------------------!
write(*,*) "Hello30"

if (associated(input%groundstate%spin)) then
write(*,*) "Hello28"
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
    if (xctype(1)==23)  then
       call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
    &          gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
    &          ec=ecir,exsr=exsr,vxsrup=vxsrup,vxsrdn=vxsrdn,v2xsrup=v2xsrup,v2xsrdn=v2xsrdn,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
       !call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, ex=exir,&
       !      ec=ecir,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
       call gv2xir_spin(gup,vxsrup,v2xsrup) 
       call gv2xir_spin(gdn,vxsrdn,v2xsrdn) 
       !call gv2xir(grho,vxsrup,v2xsrup) 
       !call gv2xir(grho,vxsrdn,v2xsrdn) 
       !CECILIA: check, but it should be fine
       vxup=vxsrup
       vxdn=vxsrdn
       exir=exsr
    elseif (xctype(1)==408)  then
       call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
    &          gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
    &          ec=ecir,exsr=exsr,vxsrup=vxsrup,vxsrdn=vxsrdn,v2xsrup=v2xsrup,v2xsrdn=v2xsrdn,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
       !call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, ex=exir,&
       !      ec=ecir,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
       call gv2xir_spin(gup,vxsrup,v2xsrup) 
       call gv2xir_spin(gdn,vxsrdn,v2xsrdn) 
       !call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, ex=exir,&
       !      ec=ecir,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
       !call gv2xir(grho,vxsr,v2xsr) 
    else
    call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
    &          gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
    &          ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
       !call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, &
       !      ex=exir,ec=ecir,vx=vx,vc=vc)
    endif
    !call xcifc(xctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup, &
    !&          gdn=gdn,g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir, &
    !&          ec=ecir,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
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
      elseif (xctype(1)==408) then !HSE
        t1=vcup(i)+vxup(i)-ex_coef*vxsrup(i)
        t2=vcdn(i)+vxdn(i)-ex_coef*vxsrdn(i)
        !vxc(1:n)=vc(1:n)+vx(1:n)-ex_coef*vxsr(1:n)     
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
      elseif (xctype(1)==408) then !HSE
        t1=vcup(i)+vxup(i)-ex_coef*vxsrup(i)
        t2=vcdn(i)+vxdn(i)-ex_coef*vxsrdn(i)
        !vxc(1:n)=vc(1:n)+vx(1:n)-ex_coef*vxsr(1:n)     
      else
        t1 = (1.d0-ex_coef)*vxup(ir)+ec_coef*vcup(ir)
        t2 = (1.d0-ex_coef)*vxdn(ir)+ec_coef*vcdn(ir)
      end if
      vxcir(ir) = 0.5d0*(t1+t2)
      bxcir(ir,1) = 0.5d0*(t1-t2)
    end do
  end if
  if (xctype(1).ne.100) then
     exir(:) = (1.d0-ex_coef)*exir(:)
  elseif (xctype(1)==408) then !HSE
     exir(:) = exir(:)-ex_coef*exsr(:)
  endif
  ecir(:) = ec_coef*ecir(:)
  
else

  !--------------------------!
  !     spin-unpolarised     !
  !--------------------------!
  if (xcgrad.le.0) then
    call xcifc(xctype,n=ngrtot,rho=rhoir,ex=exir,ec=ecir,vx=vx,vc=vc)
  else if (xcgrad.eq.1) then
    call ggair_1(grho,g2rho,g3rho)
    if (xctype(1)==23)  then
       call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, ex=exir,&
             ec=ecir,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
       call gv2xir(grho,vxsr,v2xsr) 
       vx=vxsr
       exir=exsr
    elseif (xctype(1)==408)  then
       call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, ex=exir,&
             ec=ecir,exsr=exsr,vx=vx,vc=vc,vxsr=vxsr,v2xsr=v2xsr)
       call gv2xir(grho,vxsr,v2xsr) 
    else
       call xcifc(xctype,n=ngrtot,rho=rhoir,grho=grho,g2rho=g2rho,g3rho=g3rho, &
             ex=exir,ec=ecir,vx=vx,vc=vc)
    endif
  else if (xcgrad.eq.2) then
    call ggair_2a(g2rho,gvrho,grho2)
    call xcifc(xctype,n=ngrtot,rho=rhoir,grho2=grho2,ex=exir,ec=ecir,vx=vx, &
    &          vc=vc,dxdg2=dxdg2,dcdg2=dcdg2)
    call ggair_2b(g2rho,gvrho,vx,vc,dxdg2,dcdg2)
  end if
  if ((xctype(1).ne.100) .and. (xctype(1).ne.408)) then
    vxcir(1:ngrtot) = (1.d0-ex_coef)*vx(1:ngrtot)+ec_coef*vc(1:ngrtot)
    exir(:) = (1.d0-ex_coef)*exir(:)
  elseif (xctype(1)==408) then !HSE
     vxcir(1:ngrtot)=vc(1:ngrtot)+vx(1:ngrtot)-ex_coef*vxsr(1:ngrtot)     
     exir(:) = exir(:)-ex_coef*exsr(:)
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
    if (xctype(1)==23 .or. xctype(1)==408) deallocate(exsr,vxsr,v2xsr)
  else if (xcgrad.eq.2) then
    deallocate(g2rho,gvrho,grho2)
    deallocate(dxdg2,dcdg2)
  end if
end if

return
end subroutine
!EOC

