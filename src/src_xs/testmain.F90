
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine testmain
  use modmain
  use m_getevecfvr
  use m_getevecsvr
  use m_getevalfvr
  use m_getevalsvr
  use m_getoccsvr
  implicit none
  integer, parameter :: ik=11
  integer, parameter :: istfvi=3, istfvf=7, nstfvr=istfvf-istfvi+1
  integer, parameter :: istsvi=3, istsvf=7, nstsvr=istsvf-istsvi+1
  complex(8), allocatable :: ev(:,:,:), evr(:,:,:)
  complex(8), allocatable :: evs(:,:), evsr(:,:)
  real(8), allocatable :: eva(:,:), evar(:,:)
  real(8), allocatable :: evas(:), evasr(:)
  real(8), allocatable :: occs(:), occsr(:)

  ! initialize
  call init0
  call init1

  ! first variational eigenvectors
  write(*,'(a,4i6)') 'istfvi,istfvf,nstfvr,ik',istfvi,istfvf,nstfvr,ik
  allocate(ev(nmatmax,nstfv,nspnfv),evr(nmatmax,nstfvr,nspnfv))
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),ev)
  call getevecfvr(istfvi,istfvf,vkl(:,ik),vgkl(:,:,ik,1),evr)
  write(*,'("evecfv should be zero: ",2g18.10)') ev(1,istfvi,1)-evr(1,1,1)
  write(*,'("evecfv should be zero: ",2g18.10)') ev(1,istfvf,1)-evr(1,nstfvr,1)
  write(*,'("evecfv",4g18.10)') ev(1,istfvi,1),evr(1,1,1)
  write(*,'("evecfv",4g18.10)') ev(1,istfvf,1),evr(1,nstfvr,1)
  deallocate(ev,evr)

  ! second variational eigenvectors
  write(*,'(a,4i6)') 'istsvi,istsvf,nstsvr',istsvi,istsvf,nstsvr
  allocate(evs(nstsv,nstsv),evsr(nstsvr,nstsvr))
  call getevecsv(vkl(1,ik),evs)
  call getevecsvr(istsvi,istsvf,vkl(:,ik),evsr)
  write(*,'("evecsv should be zero: ",2g18.10)') evs(istsvi,istsvi)-evsr(1,1)
  write(*,'("evecsv should be zero: ",2g18.10)') evs(istsvf,istsvi)- &
       evsr(nstsvr,1)
  write(*,'("evecsv",4g18.10)') evs(istsvi,istsvi)-evsr(1,1)
  write(*,'("evecsv",4g18.10)')evs(istsvf,istsvi)-evsr(nstsvr,1)
  deallocate(evs,evsr)

  ! first variational eigenvalues
  write(*,'(a,4i6)') 'istfvi,istfvf,nstfvr',istfvi,istfvf,nstfvr
  allocate(eva(nstfv,nspnfv),evar(nstfvr,nspnfv))
  call getevalfv(vkl(1,ik),eva)
  call getevalfvr(istfvi,istfvf,vkl(:,ik),evar)
  write(*,'("evalfv should be zero: ",g18.10)') eva(istfvi,1)-evar(1,1)
  write(*,'("evalfv should be zero: ",g18.10)') eva(istfvf,1)-evar(nstfvr,1)
  write(*,'("evalfv",2g18.10)') eva(istfvi,1),evar(1,1)
  write(*,'("evalfv",2g18.10)') eva(istfvf,1),evar(nstfvr,1)
  deallocate(eva,evar)

  ! second variational eigenvalues
  write(*,'(a,4i6)') 'istsvi,istsvf,nstsvr',istsvi,istsvf,nstsvr
  allocate(evas(nstsv),evasr(nstsvr))
  call getevalsv(vkl(1,ik),evas)
  call getevalsvr(istsvi,istsvf,vkl(:,ik),evasr)
  write(*,'("evalsv should be zero: ",g18.10)') evas(istsvi)-evasr(1)
  write(*,'("evalsv should be zero: ",g18.10)') evas(istsvf)-evasr(nstsvr)
  write(*,'("evalsv",2g18.10)') evas(istsvi),evasr(1)
  write(*,'("evalsv",2g18.10)') evas(istsvf),evasr(nstsvr)
  deallocate(evas,evasr)

  ! second variational occupation numbers
  write(*,'(a,4i6)') 'istsvi,istsvf,nstsvr',istsvi,istsvf,nstsvr
  allocate(occs(nstsv),occsr(nstsvr))
  call getoccsv(vkl(1,ik),occs)
  call getoccsvr(istsvi,istsvf,vkl(:,ik),occsr)
  write(*,'("occsv should be zero: ",g18.10)') occs(istsvi)-occsr(1)
  write(*,'("occsv should be zero: ",g18.10)') occs(istsvf)-occsr(nstsvr)
  write(*,'("occsv",2g18.10)') occs(istsvi),occsr(1)
  write(*,'("occsv",2g18.10)') occs(istsvf),occsr(nstsvr)
  deallocate(occs,occsr)

end subroutine testmain

!!$subroutine testmain
!!$  implicit none
!!$
!!$
!!$  integer,allocatable :: a1(:,:),a2(:,:), ap1(:),ap2(:), b(:,:,:,:),bp(:,:)
!!$  integer,allocatable :: bb(:,:,:,:)
!!$  integer :: i1,i2,i3,i4,j1,j2,n1,n2,n3,n4,m1,m2
!!$
!!$  n1=2; n2=3; n3=4; n4=5
!!$  m1=n1*n2; m2=n3*n4
!!$
!!$  allocate(a1(n1,n2),a2(n3,n4), ap1(m1),ap2(m2), b(n1,n2,n3,n4),bp(m1,m2))
!!$  allocate(bb(n1,n2,n3,n4))
!!$
!!$  ! assign values to arrays
!!$  do i2=1,n2
!!$     do i1=1,n1
!!$        a1(i1,i2)=5*i1+i2+1
!!$write(*,*) 'a1',i1,i2,a1(i1,i2)
!!$     end do
!!$  end do
!!$
!!$  do i4=1,n4
!!$     do i3=1,n3
!!$        a2(i3,i4)=5*i3**2+i4*3+7
!!$write(*,*) 'a2',i3,i4,a2(i3,i4)
!!$     end do
!!$  end do
!!$
!!$  
!!$  ! element-by-element
!!$  j2=0
!!$  do i4=1,n4
!!$     do i3=1,n3
!!$        j2=j2+1
!!$        j1=0
!!$        do i2=1,n2
!!$           do i1=1,n1
!!$              j1=j1+1
!!$              
!!$              b(i1,i2,i3,i4)=a1(i1,i2)*a2(i3,i4)
!!$
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$
!!$  ! pack procedure
!!$  ap1=reshape(a1,(/m1/))
!!$  ap2=reshape(a2,(/m2/))
!!$ 
!!$  do j2=1,m2
!!$     do j1=1,m1
!!$
!!$        bp(j1,j2)=ap1(j1)*ap2(j2)
!!$
!!$     end do
!!$  end do
!!$
!!$  bb=reshape(bp,(/n1,n2,n3,n4/))
!!$
!!$  do i4=1,n4
!!$     do i3=1,n3
!!$        do i2=1,n2
!!$           do i1=1,n1
!!$              
!!$              write(770077,'(4i4,4x,3i6)') i1,i2,i3,i4, &
!!$                   b(i1,i2,i3,i4),bb(i1,i2,i3,i4),b(i1,i2,i3,i4)-bb(i1,i2,i3,i4)
!!$
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$end subroutine testmain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$subroutine testmain
!!$  implicit none
!!$
!!$  logical :: ma(2,3)
!!$  integer :: f(2,3),p(6),u(2,3),j
!!$
!!$  f(1,:)=(/1,2,3/)
!!$  f(2,:)=(/7,8,9/)
!!$  p=0
!!$  p=pack(f,.true.)
!!$  ma=reshape((/(.true.,j=1,6)/),(/2,3/))
!!$  u=unpack(p,ma,0)
!!$  
!!$  write(*,*) 'f'
!!$  write(*,*) f
!!$  write(*,*) 'shape(f)'
!!$  write(*,*) shape(f)
!!$  write(*,*) 'p = packed f'
!!$  write(*,*) pack(f,.true.)
!!$  write(*,*) 'shape(p)'
!!$  write(*,*) shape(pack(f,.true.))
!!$  write(*,*) 'unpacked p'
!!$  write(*,*) unpack(p,reshape((/(.true.,j=1,6)/),(/2,3/)),0)
!!$
!!$end subroutine testmain



!!$subroutine testmain
!!$  use modmain
!!$  use m_ftfun
!!$  implicit none
!!$
!!$
!!$
!!$
!!$  complex(8), allocatable :: fmt(:,:,:)
!!$  complex(8), allocatable :: fir(:)
!!$  complex(8),allocatable :: gft(:)
!!$  integer :: ir
!!$
!!$  call init0
!!$
!!$  allocate(fmt(lmmaxvr,nrmtmax,natmtot))
!!$  allocate(fir(ngrtot),gft(ngvec))
!!$
!!$  fmt(:,:,:)=0.d0
!!$  fmt(1,:,:)=1.d0/y00
!!$  fir(:)=1.d0
!!$  call ftfun(ngvec,tir=.true.,tmt=.true.,gir=fir,gmt=fmt,ftg=gft)
!!$
!!$  do ir=1,ngvec
!!$     write(2000,'(3g18.10)') gft(ir), abs(gft(ir))
!!$  end do
!!$
!!$  write(*,*) 'normalization:',sum(abs(gft)**2)
!!$
!!$  deallocate(fir,fmt,gft)
!!$
!!$
!!$end subroutine

