
! This routine is based on code written by K. Burke.

subroutine x_pbe(kappa,mu,rho,s,u,v,ex,vx)
implicit none
! arguments
real(8), intent(in) :: kappa
real(8), intent(in) :: mu
real(8), intent(in) :: rho
real(8), intent(in) :: s
real(8), intent(in) :: u
real(8), intent(in) :: v
real(8), intent(out) :: ex
real(8), intent(out) :: vx
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: ax=-0.7385587663820224058d0
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd4=4.d0/3.d0
real(8) ul,exu,s2,p0
real(8) fxpbe,fs,fss
ul=mu/kappa
! LDA exchange energy density
exu=ax*rho**thrd
! PBE enhancement factor
s2=s**2
p0=1.d0+ul*s2
fxpbe=1.d0+kappa-kappa/p0
ex=exu*fxpbe
fs=2.d0*kappa*ul/(p0*p0)
fss=-4.d0*ul*s*fs/p0
! exchange potential
vx=exu*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
return
end subroutine

