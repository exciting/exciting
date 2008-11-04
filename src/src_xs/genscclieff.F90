
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscclieff()
  use modmain
  use modxs
  use invert
  implicit none
  ! local variables
  integer, parameter :: nleb=5810 !5810,5294
  integer :: j,itp,lm,n,ntpsph
  real(8) :: r,s,sq,m33(3,3),rnd(2)
  real(8), allocatable :: tp(:,:),spc(:,:),spcll(:,:),x(:),y(:),z(:),w(:)
  complex(8), allocatable :: ei00(:),ei00lm(:),ei00lma(:),ylm(:),zylm(:,:)
  complex(8), allocatable :: zbsht(:,:), zfsht(:,:)

  ntpsph=(lmaxspi+1)**2
  ntpsph=nleb
  allocate(ei00(ntpsph),ei00lm(ntpsph),ei00lma(lmmaxdielt))
  allocate(ylm(lmmaxdielt),zylm(ntpsph,lmmaxdielt))
  allocate(tp(2,ntpsph),spc(3,ntpsph),spcll(3,nleb))
  allocate(x(ntpsph),y(ntpsph),z(ntpsph),w(ntpsph))

  w(:)=1.d0/ntpsph
  if (lmmaxdielt.gt.ntpsph) then
     write(*,*)
     write(*,'("Error(): lmmaxdielt.lt.ntpsph: ",2i6)') lmmaxdielt,ntpsph
     write(*,*)
     stop
  end if

  ! generate Lebedev Laikov grid
  n=0
  x(:)=0.d0
  y(:)=0.d0
  z(:)=0.d0
  call ld5810(x,y,z,w,n)
  spcll(1,:)=x
  spcll(2,:)=y
  spcll(3,:)=z
 
  ! marsaglia's method
  itp=0
  do j=1,4*ntpsph
     call random_number(rnd)
     rnd=2.d0*rnd-1.d0
     s=sum(rnd**2)
     if (s.le.1.d0) then
       itp=itp+1
       sq=sqrt(1-s)
!       spcll(1,itp)=2.d0*rnd(1)*sq
!       spcll(2,itp)=2.d0*rnd(2)*sq
!       spcll(3,itp)=1.d0-2.d0*s
     end if
     if (itp.eq.ntpsph) exit
  end do


  write(*,*) 'ntpsph',ntpsph
  write(*,*) 'number of points:',n
  write(*,*) 'sum of weights',sum(w)
  write(601,*) x
  write(602,*) y
  write(603,*) z
  write(604,*) w
  write(600,'(i6,4f14.8)') (itp,x(itp),y(itp),z(itp),w(itp),itp=1,nleb)

  ! preset dielectric tensor for testing
  dielten(1,:)=(/ 2.91911039, 0.00000000, 3.49765354 /)
  dielten(2,:)=(/ 0.00000000, 2.79383654, 0.00000000 /)
  dielten(3,:)=(/ 3.49765354, 0.00000000, 102.25001110 /)
  dielten(1,:)=dielten(1,:)+zi*(/ 0.00000000, 0.00000000, 0.00000579 /)
  dielten(2,:)=dielten(2,:)+zi*(/ 0.00000000, 0.00000000, 0.00000000 /)
  dielten(3,:)=dielten(3,:)+zi*(/ -0.00000579,0.00000000, 0.00000000 /)


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  call random_number(m33)
  dielten(:,:)=dielten(:,:)+m33(:,:)*3.d0
  call random_number(m33)
  dielten(:,:)=dielten(:,:)+zi*m33(:,:)*3.d0
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! generate spherical covering set (angles and coordinates)
  call sphcover(ntpsph,tp)
  spc(1,:)=sin(tp(1,:))*cos(tp(2,:))
  spc(2,:)=sin(tp(1,:))*sin(tp(2,:))
  spc(3,:)=cos(tp(1,:))
  
  spc=spcll !********
  write(*,*) 'norm of vectors:',sum(spc**2)
  
  ! generate spherical harmonics on large covering set
  do itp=1,ntpsph
     call sphcrd(spc(:,itp),r,tp(:,itp))
     call genylm(lmaxdielt,tp(:,itp),ylm)
     zylm(itp,:)=ylm(:)
  end do

  ! calculate function on covering set
  do itp=1,ntpsph
     ! head which is the 1/(p*L*p) factor
     ei00(itp)=1.d0/dot_product(spc(:,itp),matmul(dielten,spc(:,itp)))
  end do

  ! calculate lm-expansion coefficients
  do lm=1,lmmaxdielt
     ei00lma(lm)=fourpi*dot_product(zylm(:,lm),ei00*w)
  end do

  write(*,*) 'dielten',dielten
  write(10000+lmaxspi,'(i6,3f14.8)') (itp,spc(:,itp),itp=1,ntpsph)
  write(20000+lmaxspi,'(i6,2f14.8)') (itp,ei00(itp),itp=1,ntpsph)
  write(40000+lmaxspi,'(i6,2f14.8)') (lm,ei00lma(lm),lm=1,lmmaxdielt)

  deallocate(ei00,ei00lm,ei00lma,ylm,zylm,tp,spc,spcll,x,y,z,w)

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
