
! Copyright (C) 2006 Zhigang Wu and R. E. Cohen.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine xc_wc06(n,rho,grho,g2rho,g3rho,ex,ec,vx,vc)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rho(n)
real(8), intent(in) :: grho(n)
real(8), intent(in) :: g2rho(n)
real(8), intent(in) :: g3rho(n)
real(8), intent(out) :: ex(n)
real(8), intent(out) :: ec(n)
real(8), intent(out) :: vx(n)
real(8), intent(out) :: vc(n)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd2=2.d0/3.d0
! default PBE beta
real(8), parameter :: beta=0.06672455060314922d0
! maximum allowed |grad rho|
real(8), parameter :: gmax=1.d6
! maximum allowed grad^2 rho
real(8), parameter :: g2max=1.d12
! maximum allowed (grad rho).(grad |grad rho|)
real(8), parameter :: g3max=1.d14
real(8) r,grho_,g2rho_,g3rho_
real(8) kf,s,u,v,rs,z,g
real(8) ks,ksg,t,uu,vv,ww
do i=1,n
  r=rho(i)
  if (r.gt.1.d-12) then
    grho_=grho(i)
    g2rho_=g2rho(i)
    g3rho_=g3rho(i)
! check gradients are within range
    if (grho_.gt.gmax) grho_=gmax
    if (abs(g2rho_).gt.g2max) g2rho_=sign(g2max,g2rho_)
    if (abs(g3rho_).gt.g3max) g3rho_=sign(g3max,g3rho_)
    kf=(r*3.d0*pi**2)**thrd
    s=grho_/(2.d0*kf*r)
    u=g3rho_/((r**2)*(2.d0*kf)**3)
    v=g2rho_/(r*(2.d0*kf)**2)
! Wu-Cohen exchange
    call x_wc06(r,s,u,v,ex(i),vx(i))
! Perdew-Burke-Ernzerhof correlation
    rs=(3.d0/(4.d0*pi*r))**thrd
    z=0.d0
    g=1.d0
    ks=sqrt(4.d0*kf/pi)
    ksg=2.d0*ks*g
    t=grho_/(ksg*r)
    uu=g3rho_/((r**2)*ksg**3)
    vv=g2rho_/(r*ksg**2)
    ww=0.d0
    call c_pbe(beta,rs,z,t,uu,vv,ww,ec(i),vc(i),vc(i))
  else
    ex(i)=0.d0
    ec(i)=0.d0
    vx(i)=0.d0
    vc(i)=0.d0
  end if
end do
return
end subroutine

