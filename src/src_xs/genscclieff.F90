
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscclieff()
  use modmain
  use modxs
  use invert
  implicit none
  ! local variables
  integer, parameter :: nsphcov=5810
  integer :: itp,lm,ntpsph
  real(8) :: t1,t2,r,qsz,clwt,clwt2, ts0,ts1
  complex(8) :: z00,z01,zt1,zt2
  real(8), allocatable :: plat(:,:),p(:),tp(:,:),spc(:,:),w(:)
  complex(8), allocatable :: m00lm(:),mx0lm(:),mxxlm(:)
  complex(8), allocatable :: ylm(:),zylm(:,:)
  complex(8), allocatable :: ei00(:),ei00lm(:),ei00lma(:)
  ! *** values for PA ***
  call preset_dielten
  zt1=1.d0/((dielten(1,1)+dielten(2,2)+dielten(3,3))/3.d0)
  zt2=(1.d0/dielten(1,1)+1.d0/dielten(2,2)+1.d0/dielten(3,3))/3.d0
  ! Wigner-Seitz radius and spherical approximation to 1/q^2 average
  qsz=(6*pi**2/(omega*product(ngridq)))**(1.d0/3.d0)
  clwt=2*qsz*omega*product(ngridq)/pi
!!$  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$  call init2
!!$  call genwiq2
!!$  call writewiq2
!!$  clwt2=wiq2(1)*fourpi*product(ngridq)*omega/(twopi**3)
!!$  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! number of points on sphere
  if (tleblaik) then
     ntpsph=nleblaik
  else
     ntpsph=nsphcov
  end if
  if (lmmaxdielt.gt.ntpsph) then
     write(*,*)
     write(*,'("Error(genscclieff): lmmaxdielt.gt.ntpsph: ",2i6)') lmmaxdielt, &
          ntpsph
     write(*,*)
     stop
  end if
  allocate(plat(3,ntpsph),p(ntpsph))
  allocate(m00lm(lmmaxdielt),mx0lm(lmmaxdielt),mxxlm(lmmaxdielt))
  allocate(ei00(ntpsph),ei00lm(ntpsph),ei00lma(lmmaxdielt))
  allocate(ylm(lmmaxdielt),zylm(ntpsph,lmmaxdielt))
  allocate(tp(2,ntpsph),spc(3,ntpsph))
  allocate(w(ntpsph))
  if (tleblaik) then
     ! generate Lebedev Laikov grid
     call leblaik(ntpsph,spc,w)
     ! generate tetha and phi angles
     do itp=1,ntpsph
        call sphcrd(spc(:,itp),r,tp(:,itp))
     end do
  else
     ! distribution is assumed to be uniform
     w(:)=1.d0/ntpsph
     ! generate spherical covering set (angles and coordinates)
     call sphcover(ntpsph,tp)
     spc(1,:)=sin(tp(1,:))*cos(tp(2,:))
     spc(2,:)=sin(tp(1,:))*sin(tp(2,:))
     spc(3,:)=cos(tp(1,:))
  end if    
  ! generate spherical harmonics on covering set
  do itp=1,ntpsph
     call genylm(lmaxdielt,tp(:,itp),ylm)
     zylm(itp,:)=ylm(:)
  end do

  ! unit vectors of spherical covering set in lattice coordinates
  plat=matmul(binv,spc)
  ! distances to subcell cell boundaries in reciprocal space
  do itp=1,ntpsph
     p(itp:)=1.d0/(2.d0*maxval(abs(ngridq(:)*plat(:,itp)),1))
  end do

  ! calculate function on covering set
  do itp=1,ntpsph
     ! head which is the 1/(p*L*p) factor
     ei00(itp)=1.d0/dot_product(spc(:,itp),matmul(dielten,spc(:,itp)))
  end do

  ! calculate lm-expansion coefficients
  do lm=1,lmmaxdielt
     ei00lma(lm)=fourpi*dot_product(zylm(:,lm),ei00*w)
     m00lm(lm)=fourpi*dot_product(zylm(:,lm),p*w)
  end do

  ! scaling factor
  t1=(omega/(twopi)**3)*product(ngridq)*fourpi
  ! spherical approximation
  z00=t1*dot_product(m00lm,ei00lma)
  ! subcell average
  z01=t1*conjg(m00lm(1))*ei00lma(1)

  t2=product(ngridq)*omega
  write(*,*) 'analytic average (av of scr):', zt1*clwt/t2
  write(*,*) 'analytic average (av of inv.scr):',zt2*clwt/t2
  write(*,*) 'spherical average',z00/t2
  write(*,*) 'spherical average(1)',z01/t2
  write(*,*) 'potcl eff',clwt
  write(*,*) 'potcl eff2',clwt2
  write(*,*) 'qsz',qsz
  write(*,*) 'V_gamma',fourpi*qsz**3/3.d0
  write(*,*) 'z00',z00
  write(*,*) 'z01',z01
  write(*,*) 'sph.appr.',1.d0/(2.d0*pi)**3*fourpi**2 *qsz *t2
  write(*,*) 'm00lm(1)',m00lm(1)
  write(*,*) 'ei00lma(1)',ei00lma(1)
  write(*,*) 'm00lm(1)*ei00lma(1)',ei00lma(1)*m00lm(1)
  write(*,*) 'dielten',dielten
  write(10000+lmaxdielt,'(i6,3f14.8)') (itp,spc(:,itp),itp=1,ntpsph)
  write(20000+lmaxdielt,'(i6,2f14.8)') (itp,ei00(itp),itp=1,ntpsph)
  write(40000+lmaxdielt,'(i6,2f14.8)') (lm,ei00lma(lm),lm=1,lmmaxdielt)
  write(45000+lmaxdielt,'(i6,2f14.8)') (lm,m00lm(lm),lm=1,lmmaxdielt)
  write(50000+lmaxdielt,'(2i6,3g18.10)') lmaxdielt,lmmaxdielt, &
       dble(z00)/t2,dble(zt1*clwt/t2),dble(zt2*clwt/t2)

  deallocate(ei00,ei00lm,ei00lma,m00lm,mx0lm,mxxlm,ylm,zylm,tp,spc,w,plat,p)

end subroutine genscclieff





subroutine angavdm(n,l,s,bi,av)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: n
  complex(8), intent(in) :: bi(n-1,n-1),s(n-1,3),l(3,3)
  complex(8), intent(out) :: av(n,n)
  ! local variables
  integer :: i,j,lm
  real(8) :: t1
  real(8), allocatable :: plat(:,:),p(:)
  complex(8), allocatable :: m00lm(:),m0xlm(:)
  complex(8), allocatable :: ei00(:),ei0x(:,:),ei00lm(:),ei0xlm(:,:)
  allocate(plat(3,lmmaxvr),p(lmmaxvr))
  allocate(m00lm(lmmaxvr),m0xlm(lmmaxvr))
  allocate(ei00(lmmaxvr),ei0x(lmmaxvr,n-1))
  allocate(ei00lm(lmmaxvr),ei0xlm(lmmaxvr,n-1))
  ! unit vectors of spherical covering set in lattice coordinates
  plat=matmul(binv,sphcov)
  ! distances to unit cell boundary in reciprocal space
  p(:)=1.d0/(2.d0*maxval(abs(plat),1))
  do lm=1,lmmaxvr
     ! head which is the 1/(p*L*p) factor
     ei00(lm)=1.d0/dot_product(sphcov(:,lm),matmul(l,sphcov(:,lm)))
     ! wings: -p*s * ei00
     ei0x(lm,:)=-ei00(lm)*matmul(s,sphcov(:,lm))
  end do
! spherical approximation
!***p(:)=((3.d0*(twopi**3/omega))/(fourpi*product(ngridq)))**(1.d0/3.d0)
  ! forward transform shape dependent factor for head to spherical harmonics
  m00lm(:)=matmul(zfshtvr,p)
  ! forward transform shape dependent factor for wings to spherical harmonics
  m0xlm(:)=matmul(zfshtvr,(1.d0/2.d0)*p**2)
  !forward transform to spherical harmonics
  ei00lm=matmul(zfshtvr,ei00)
  ei0xlm=matmul(zfshtvr,ei0x)
  ! angular average
  t1=(omega/(twopi)**3)/product(ngridq)
  av(1,1)=t1*dot_product(m00lm,ei00lm)
  av(1,2:)=t1*matmul(m0xlm,ei0xlm)
  ! body without angular average
  av(2:,2:)=bi(:,:)
  ! set lower triangle
  do i=1,n
     do j=i+1,n
        av(j,i)=conjg(av(i,j))
     end do
  end do 
  deallocate(plat,p,m00lm,m0xlm,ei00,ei0x,ei00lm,ei0xlm)
  
  write(7788,'(2i8,3g18.10)') ((i,j,av(i,j),abs(av(i,j))**2,j=1,n),i=1,n)
   
end subroutine angavdm


subroutine preset_dielten
  use modmain
  use modxs
  implicit none
  real(8) :: r(3,3)
  ! preset dielectric tensor for testing
  dielten(:,:)=zzero
!!$!  (values are for trans-polyacetylene) from 2x2x16 k-point grid
!!$  dielten(1,:)=(/ 2.91911039, 0.00000000, 3.49765354 /)
!!$  dielten(2,:)=(/ 0.00000000, 2.79383654, 0.00000000 /)
!!$  dielten(3,:)=(/ 3.49765354, 0.00000000, 102.25001110 /)
!!$  dielten(1,:)=dielten(1,:)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
!!$  dielten(2,:)=dielten(2,:)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
!!$  dielten(3,:)=dielten(3,:)+zi*(/ -0.00000579,0.00000000, 0.00000000 /)
!  (values are for trans-polyacetylene), 4x4x32 k-point grid
  dielten(1,:)=(/ 2.91911039, 0.00000000, 0.49765354 /)
  dielten(2,:)=(/ 0.00000000, 2.79383654, 0.00000000 /)
  dielten(3,:)=(/ 0.49765354, 0.00000000, 52.725001110 /)
  dielten(1,:)=dielten(1,:)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
  dielten(2,:)=dielten(2,:)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
  dielten(3,:)=dielten(3,:)+zi*(/ -0.00000579,0.00000000, 0.00000000 /)
!!$! testing only
!!$  dielten(1,:)=(/ 3.0, 0.0, 0.0 /)
!!$  dielten(2,:)=(/ 0.0, 3.0, 0.0 /)
!!$  dielten(3,:)=(/ 0.0, 0.0, 3.0 /)
!!$  dielten(1,:)=dielten(1,:)+zi*(/ 0.0, 0.0, 0.0 /)
!!$  dielten(2,:)=dielten(2,:)+zi*(/ 0.0, 0.0, 0.0 /)
!!$  dielten(3,:)=dielten(3,:)+zi*(/ 0.0, 0.0, 0.0 /)
!!$  call random_number(r)
!!$  dielten(:,:)=dielten(:,:)+r(:,:)*1.d0
!!$  call random_number(r)
!!$  dielten(:,:)=dielten(:,:)+zi*r(:,:)*1.d0
end subroutine preset_dielten
