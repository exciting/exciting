
! Copyright (C) 2006 Zhigang Wu and R. E. Cohen.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine x_wc06(rho,s,u,v,ex,vx)
implicit none
! arguments
real(8), intent(in) :: rho
real(8), intent(in) :: s
real(8), intent(in) :: u
real(8), intent(in) :: v
real(8), intent(out) :: ex
real(8), intent(out) :: vx
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: ax=-0.7385587663820224059d0
real(8), parameter :: mu=0.2195149727645171d0
real(8), parameter :: kappa=0.804d0
real(8), parameter :: b=10.d0/81.d0
real(8), parameter :: c=0.00793746933516d0
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd4=4.d0/3.d0
real(8) dmu,exu
real(8) s2,s4,es2,x,p0,fxwc
real(8) fs,fss,t0,t1,t2,t3
! lda exchange energy density
exu=ax*rho**thrd
s2=s**2
s4=s2**2
es2=exp(-s2)
t0=1.d0+c*s4
dmu=mu-b
x=b*s2+dmu*s2*es2+log(t0)
p0=1.d0+x/kappa
! WC enhancement factor
fxwc=1.d0+kappa-kappa/p0
! exchange energy density
ex=exu*fxwc
t1=b+dmu*(1.d0-s2)*es2+2.d0*c*s2/t0
t2=dmu*s*(s2-2.d0)*es2+2.d0*c/t0-4.d0*(c**2)*s4/(t0**2)
t3=1.d0/(p0**2)
fs=2.d0*t1*t3
fss=t3*(4.d0*t2-8.d0*s*(t1**2)/(kappa*p0))
! exchange potential
vx=exu*(thrd4*fxwc-(u-thrd4*s2*s)*fss-v*fs)
return
end

