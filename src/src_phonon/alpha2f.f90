
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine alpha2f
use modmain
implicit none
! local variables
integer n,ik,iq,i,j
integer i1,i2,i3,iw
integer lwork,info
real(8) wmin,wmax,wd,dw,wlog
real(8) v(3),lambda,tc,t1
! allocatable arrays
real(8), allocatable :: wq(:,:)
real(8), allocatable :: wp(:)
real(8), allocatable :: gq(:,:)
real(8), allocatable :: a2fp(:)
real(8), allocatable :: w(:)
real(8), allocatable :: a2f(:)
real(8), allocatable :: f(:),g(:),cf(:,:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: dynp(:,:)
complex(8), allocatable :: dynr(:,:,:)
complex(8), allocatable :: ev(:,:),b(:,:)
complex(8), allocatable :: a2fmq(:,:,:)
complex(8), allocatable :: a2fmr(:,:,:)
complex(8), allocatable :: a2fmp(:,:)
complex(8), allocatable :: work(:)
! initialise universal variables
call init0
call init1
call init2
n=3*natmtot
allocate(wq(n,nqpt))
allocate(wp(n))
allocate(gq(n,nqpt))
allocate(a2fp(n))
allocate(w(nwdos))
allocate(a2f(nwdos))
allocate(f(nwdos),g(nwdos),cf(3,nwdos))
allocate(rwork(3*n))
allocate(dynq(n,n,nqpt))
allocate(dynp(n,n))
allocate(dynr(n,n,ngridq(1)*ngridq(2)*ngridq(3)))
allocate(ev(n,n),b(n,n))
allocate(a2fmq(n,n,nqpt))
allocate(a2fmr(n,n,ngridq(1)*ngridq(2)*ngridq(3)))
allocate(a2fmp(n,n))
lwork=2*n
allocate(work(lwork))
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! compute the density of states at the Fermi energy
call occupy
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! read in the phonon linewidths for each q-point
call readgamma(gq)
! loop over phonon q-points
do iq=1,nqpt
! diagonalise the dynamical matrix
  call dyndiag(dynq(:,:,iq),wq(:,iq),ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues are the phonon linewidths divided by the frequency
  do i=1,n
    if (wq(i,iq).gt.1.d-8) then
      t1=gq(i,iq)/wq(i,iq)
    else
      t1=0.d0
    end if
    do j=1,n
      b(i,j)=t1*conjg(ev(j,i))
    end do
  end do
  call zgemm('N','N',n,n,n,zone,ev,n,b,n,zzero,a2fmq(:,:,iq),n)
end do
! Fourier transform the matrices to real-space
call dynqtor(a2fmq,a2fmr)
! find the minimum and maximum frequencies
wmin=0.d0
wmax=0.d0
do iq=1,nqpt
  wmin=min(wmin,wq(1,iq))
  wmax=max(wmax,wq(n,iq))
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
a2f(:)=0.d0
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
! compute the alpha^2F matrix at this particular q-point
      call dynrtoq(v,a2fmr,a2fmp)
! diagonlise the alpha^2F matrix
      call zheev('N','U',n,a2fmp,n,a2fp,work,lwork,rwork,info)
      do i=1,n
        t1=(wp(i)-wmin)/dw+1.d0
        iw=nint(t1)
        if ((iw.ge.1).and.(iw.le.nwdos)) then
          a2f(iw)=a2f(iw)+a2fp(i)
        end if
      end do
    end do
  end do
end do
t1=twopi*(fermidos/2.d0)*dw*dble(ngrdos)**3
if (t1.gt.1.d-8) then
  t1=1.d0/t1
else
  t1=0.d0
end if
a2f(:)=t1*a2f(:)
! smooth Eliashberg function if required
if (nsmdos.gt.0) call fsmooth(nsmdos,nwdos,1,a2f)
! write Eliashberg function to file
open(50,file='ALPHA2F.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),a2f(iw)
end do
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" Eliashberg function written to ALPHA2F.OUT")')
! compute the total lambda
do iw=1,nwdos
  if (w(iw).gt.1.d-8) then
    f(iw)=a2f(iw)/w(iw)
  else
    f(iw)=0.d0
  end if
end do
call fderiv(-1,nwdos,w,f,g,cf)
lambda=2.d0*g(nwdos)
open(50,file='LAMBDA.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("Electron-phonon mass enhancement parameter, lambda : ",G18.10)') &
 lambda
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" Electron-phonon mass enhancement parameter, lambda, written to&
 & LAMBDA.OUT")')
! compute the logarithmic average frequency
do iw=1,nwdos
  if (w(iw).gt.1.d-8) then
    f(iw)=a2f(iw)*log(w(iw))/w(iw)
  else
    f(iw)=0.d0
  end if
end do
call fderiv(-1,nwdos,w,f,g,cf)
t1=(2.d0/lambda)*g(nwdos)
wlog=exp(t1)
! compute McMillan-Allen-Dynes superconducting critical temperature
t1=(-1.04d0*(1.d0+lambda))/(lambda-mustar-0.62d0*lambda*mustar)
tc=(wlog/(1.2d0*kboltz))*exp(t1)
open(50,file='TC_MCMILLAN.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("Logarithmic average frequency : ",G18.10)') wlog
write(50,*)
write(50,'("Coulomb pseudopotential, mu* : ",G18.10)') mustar
write(50,*)
write(50,'("McMillan-Allen-Dynes superconducting critical temperature")')
write(50,'("[Eq. 34, Phys. Rev. B 12, 905 (1975)] (kelvin) : ",G18.10)') tc
write(50,*)
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" McMillan-Allen-Dynes superconducting critical temperature, T_c,&
 & written to TC_MCMILLAN.OUT")')
write(*,*)
deallocate(wq,wp,gq,a2fp,w,a2f,f,g,cf)
deallocate(rwork,dynq,dynp,dynr,ev,b)
deallocate(a2fmq,a2fmr,a2fmp,work)
return
end subroutine

