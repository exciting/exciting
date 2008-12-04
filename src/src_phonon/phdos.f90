
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdos
use modmain
implicit none
! local variables
! number of temperature values
integer, parameter :: ntemp=10
integer n,iq,i,iw
integer i1,i2,i3
real(8) wmin,wmax,wd,dw
real(8) tmax,temp(ntemp),s(ntemp)
real(8) v(3),t1,t2
! allocatable arrays
real(8), allocatable :: wp(:)
real(8), allocatable :: w(:)
real(8), allocatable :: gw(:)
real(8), allocatable :: f(:),g(:),cf(:,:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: dynr(:,:,:)
complex(8), allocatable :: dynp(:,:)
complex(8), allocatable :: ev(:,:)
! initialise universal variables
call init0
call init2
n=3*natmtot
allocate(wp(n))
allocate(w(nwdos))
allocate(gw(nwdos))
allocate(f(nwdos),g(nwdos),cf(3,nwdos))
allocate(dynq(n,n,nqpt))
allocate(dynr(n,n,ngridq(1)*ngridq(2)*ngridq(3)))
allocate(dynp(n,n))
allocate(ev(n,n))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! find the minimum and maximum frequencies
wmin=0.d0
wmax=0.d0
do iq=1,nqpt
  call dyndiag(dynq(:,:,iq),wp,ev)
  wmin=min(wmin,wp(1))
  wmax=max(wmax,wp(n))
end do
wmax=wmax+(wmax-wmin)*0.1d0
wmin=wmin-(wmax-wmin)*0.1d0
wd=wmax-wmin
if (wd.lt.1.d-8) wd=1.d0
dw=wd/dble(nwdos)
! generate energy grid
do iw=1,nwdos
  w(iw)=dw*dble(iw-1)+wmin
end do
gw(:)=0.d0
do i1=0,ngrdos-1
  v(1)=dble(i1)/dble(ngrdos)
  do i2=0,ngrdos-1
    v(2)=dble(i2)/dble(ngrdos)
    do i3=0,ngrdos-1
      v(3)=dble(i3)/dble(ngrdos)
! compute the dynamical matrix at this particular q-point
      call dynrtoq(v,dynr,dynp)
! find the phonon frequencies
      call dyndiag(dynp,wp,ev)
      do i=1,n
        t1=(wp(i)-wmin)/dw+1.d0
        iw=nint(t1)
        if ((iw.ge.1).and.(iw.le.nwdos)) then
          gw(iw)=gw(iw)+1.d0
        end if
      end do
    end do
  end do
end do
t1=1.d0/(dw*dble(ngrdos)**3)
gw(:)=t1*gw(:)
! smooth phonon DOS if required
if (nsmdos.gt.0) call fsmooth(nsmdos,nwdos,1,gw)
! write phonon DOS to file
open(50,file='PHDOS.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),gw(iw)
end do
close(50)
write(*,*)
write(*,'("Info(phdos):")')
write(*,'(" phonon density of states written to PHDOS.OUT")')
!-------------------------------------------!
!     thermodynamic properties from DOS     !
!-------------------------------------------!
! maximum temperature
tmax=wmax/kboltz
! temperature grid
do i=1,ntemp
  temp(i)=tmax*dble(i)/dble(ntemp)
end do
open(50,file='THERMO.OUT',action='WRITE',form='FORMATTED')
! zero point energy
do iw=1,nwdos
  f(iw)=gw(iw)*w(iw)
end do
call fderiv(-1,nwdos,w,f,g,cf)
t1=0.5d0*dble(natmtot)*g(nwdos)
write(50,*)
write(50,'("Zero-point energy : ",G18.10)') t1
! vibrational energy
write(50,*)
write(50,'("Vibrational energy vs. temperature :")')
do i=1,ntemp
  do iw=1,nwdos
    t1=w(iw)/(2.d0*kboltz*temp(i))
    if (t1.gt.0.d0) then
      f(iw)=gw(iw)*w(iw)*cosh(t1)/sinh(t1)
    else
      f(iw)=0.d0
    end if
  end do
  call fderiv(-1,nwdos,w,f,g,cf)
  t1=0.5d0*dble(natmtot)*g(nwdos)
  write(50,'(2G18.10)') temp(i),t1
  s(i)=t1
end do
! free energy
write(50,*)
write(50,'("Free energy vs. temperature :")')
do i=1,ntemp
  do iw=1,nwdos
    t1=2.d0*sinh(w(iw)/(2.d0*kboltz*temp(i)))
    if (t1.gt.0.d0) then
      f(iw)=gw(iw)*log(t1)
    else
      f(iw)=0.d0
    end if
  end do
  call fderiv(-1,nwdos,w,f,g,cf)
  t1=dble(natmtot)*kboltz*temp(i)*g(nwdos)
  write(50,'(2G18.10)') temp(i),t1
! compute entropy from S = (F-E)/T
  s(i)=(s(i)-t1)/temp(i)
end do
! entropy
write(50,*)
write(50,'("Entropy vs. temperature :")')
do i=1,ntemp
  write(50,'(2G18.10)') temp(i),s(i)
end do
! heat capacity
write(50,*)
write(50,'("Heat capacity vs. temperature :")')
do i=1,ntemp
  do iw=1,nwdos
    t1=w(iw)/(kboltz*temp(i))
    t2=exp(t1)-1.d0
    if (t2.ne.0.d0) then
      f(iw)=gw(iw)*(t1**2)*(t2+1.d0)/t2**2
    else
      f(iw)=0.d0
    end if
  end do
  call fderiv(-1,nwdos,w,f,g,cf)
  t1=dble(natmtot)*kboltz*g(nwdos)
  write(50,'(2G18.10)') temp(i),t1
end do
close(50)
write(*,'(" thermodynamic properties written to THERMO.OUT")')
write(*,*)
deallocate(wp,w,gw,f,g,cf,dynq,dynr,dynp,ev)
return
end subroutine

