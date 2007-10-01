
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdos
use modmain
implicit none
! local variables
integer n,iq,i,i1,i2,i3,iw
real(8) wmin,wmax,wd,dw,h,v(3),t1
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: g(:)
complex(8), allocatable :: ev(:,:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: dynp(:,:)
complex(8), allocatable :: dynr(:,:,:)
! initialise universal variables
call init0
call init2
n=3*natmtot
allocate(w(n))
allocate(g(nwdos))
allocate(ev(n,n))
allocate(dynq(n,n,nqpt))
allocate(dynp(n,n))
allocate(dynr(n,n,ngridq(1)*ngridq(2)*ngridq(3)))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! find the minimum and maximum frequencies
wmin=1.d8
wmax=0.d0
do iq=1,nqpt
  call dyndiag(dynq(1,1,iq),w,ev)
  wmin=min(wmin,w(1))
  wmax=max(wmax,w(n))
end do
wmax=wmax+(wmax-wmin)*0.1d0
wmin=wmin-(wmax-wmin)*0.1d0
wd=wmax-wmin
dw=wd/dble(nwdos)
h=1.d0/dw
g(:)=0.d0
do i1=0,ngrdos-1
  v(1)=dble(i1)/dble(ngrdos)
  do i2=0,ngrdos-1
    v(2)=dble(i2)/dble(ngrdos)
    do i3=0,ngrdos-1
      v(3)=dble(i3)/dble(ngrdos)
      call dynrtoq(v,dynr,dynp)
      call dyndiag(dynp,w,ev)
      do i=1,n
        t1=(w(i)-wmin)/dw+1.d0
        iw=nint(t1)
        if ((iw.ge.1).and.(iw.le.nwdos)) then
          g(iw)=g(iw)+h
        end if
      end do
    end do
  end do
end do
t1=1.d0/(dble(ngrdos)**3)
g(:)=t1*g(:)
! smooth phonon DOS if required
if (nsmdos.gt.0) call fsmooth(nsmdos,nwdos,1,g)
! write phonon DOS to file
open(50,file='PHDOS.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  t1=dble(iw-1)*dw+wmin
  write(50,'(2G18.10)') t1,g(iw)
end do
close(50)
write(*,*)
write(*,'("Info(phdos):")')
write(*,'(" phonon density of states written to PHDOS.OUT")')
write(*,*)
deallocate(w,g,ev,dynq,dynp,dynr)
return
end subroutine

