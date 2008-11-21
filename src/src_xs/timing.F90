
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xstiming
  use invert
  use blaswrappers
  use m_getunit
  implicit none
  ! local variables
  integer :: un,nlo,nup,ninc,j,k,npts
  real(8) :: cpu0,cpu1,ts0,ts1,tc_mm,tc_inv,tc_diag,ts_mm,ts_inv,ts_diag
  real(8) :: f
  complex(8) :: alpha,beta
  real(8), allocatable :: v(:),mr(:,:),tc(:),x(:),c(:)
  complex(8), allocatable :: mz(:,:),mh(:,:),ma(:,:),mb(:,:)
  real(8), external :: polynom
  nlo=100
  nup=1000
  ninc=100
  alpha=(2.d0,3.d0)
  beta=(4.d0,5.d0)
  call getunit(un)
  open(un,file='TIMING.OUT',form='formatted',action='write',status='replace')
  npts=(nup-nlo)/ninc+1
  allocate(x(npts),tc(npts),c(npts))
  k=0
  do j=nlo,nup,ninc
     k=k+1
     allocate(mr(j,j),v(j),mz(j,j),mh(j,j),ma(j,j),mb(j,j))
     ! set up random number matrix
     mz(:,:)=0.d0
     call random_number(mr)
     mz=mz+mr
     call random_number(mr)
     mz=mz+(0.d0,1.d0)*mr
     ma=mz
     ! test matrix multiplication
     call timesec(ts0)
     call cpu_time(cpu0)
     call zgemm_wrap(mb,'n',ma,'n',mz,alpha,beta)
     call timesec(ts1)
     call cpu_time(cpu1)
     ts_mm=ts1-ts0
     tc_mm=cpu1-cpu0
     ! test inversion
     call timesec(ts0)
     call cpu_time(cpu0)
     call zinvert_lapack(mz,mb)
     call timesec(ts1)
     call cpu_time(cpu1)
     ts_inv=ts1-ts0
     tc_inv=cpu1-cpu0
     ! test diagonalization
     call timesec(ts0)
     call cpu_time(cpu0)
     call bsesoldiag(j,j,mz,v,mb)
     call timesec(ts1)
     call cpu_time(cpu1)
     ts_diag=ts1-ts0
     tc_diag=cpu1-cpu0
     write(un,'(i8,3g18.10)') j,tc_mm,tc_inv,tc_diag
     deallocate(mr,v,mz,mh,ma,mb)
     x(k)=j
     tc(k)=tc_diag
  end do
  close(un)
  f=polynom(0,5,x,tc,c,1000.d0)
  write(*,*) 'polynomial fit: coefficients a0,a1,a2,a3',c
  write(*,*) 'extrapolation to 1000:',f
end subroutine xstiming
