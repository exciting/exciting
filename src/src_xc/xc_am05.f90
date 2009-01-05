
! Copyright (C) 2004, 2005 Rickard Armiento
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: xc_am05
! !INTERFACE:
subroutine xc_am05(n,rho,grho,g2rho,g3rho,ex,ec,vx,vc)
! !INPUT/OUTPUT PARAMETERS:
!   n     : number of density points (in,integer)
!   rho   : charge density (in,real(n))
!   grho  : |grad rho| (in,real(n))
!   g2rho : grad^2 rho (in,real(n))
!   g3rho : (grad rho).(grad |grad rho|) (in,real(n))
!   ex    : exchange energy density (out,real(n))
!   ec    : correlation energy density (out,real(n))
!   vx    : spin-unpolarised exchange potential (out,real(n))
!   vc    : spin-unpolarised correlation potential (out,real(n))
! !DESCRIPTION:
!   Spin-unpolarised exchange-correlation potential and energy functional of
!   R. Armiento and A. E. Mattsson, {\it Phys. Rev. B} {\bf 72}, 085108 (2005).
!
! !REVISION HISTORY:
!   Created April 2005 (RAR); based on xc_pbe
!EOP
!BOC
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
! maximum allowed |grad rho|
real(8), parameter :: gmax=1.d6
! maximum allowed grad^2 rho
real(8), parameter :: g2max=1.d12
! maximum allowed (grad rho).(grad |grad rho|)
real(8), parameter :: g3max=1.d14
real(8) r,kf,s,v,u
real(8) grho_,g2rho_,g3rho_
do i=1,n
  if (rho(i).gt.1.d-12) then
    grho_=grho(i)
    g2rho_=g2rho(i)
    g3rho_=g3rho(i)
! check gradients are within range
    if (grho_.gt.gmax) grho_=gmax
    if (abs(g2rho_).gt.g2max) g2rho_=sign(g2max,g2rho_)
    if (abs(g3rho_).gt.g3max) g3rho_=sign(g3max,g3rho_)
! exchange energy density and potential
    r=rho(i)
    kf=(r*3.d0*pi**2)**(1.d0/3.d0)
    s=grho_/(2.d0*kf*r)
    v=g2rho_/(r*(2.d0*kf)**2)
    u=g3rho_/((r**2)*(2.d0*kf)**3)
    call xc_am05_point(r,s,u,v,ex(i),ec(i),vx(i),vc(i),1)
  else
    ex(i)=0.d0
    ec(i)=0.d0
    vx(i)=0.d0
    vc(i)=0.d0
  end if
end do
return
end subroutine
!EOC

!BOP
! !ROUTINE: xc_am05_point
! !INTERFACE:
subroutine xc_am05_point(rho,s,u,v,ex,ec,vx,vc,pot)
! !INPUT/OUTPUT PARAMETERS:
!   rho : electron density (in,real)
!   s   : gradient of n / (2 kF n)
!   u   : grad n * grad | grad n | / (n**2 (2 kF)**3)
!   v   : laplacian of density / (n**2 (2.d0*kf)**3)
!   ex  : exchange energy density (out,real)
!   ec  : correlation energy density (out,real)
!   vx  : spin-unpolarised exchange potential (out,real)
!   vc  : spin-unpolarised correlation potential (out,real)
! !DESCRIPTION:
!   Calculate the spin-unpolarised exchange-correlation potential and energy for
!   the Armiento-Mattsson 05 functional for a single point.
!
! !REVISION HISTORY:
!   Created April 2005 (RAR)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rho, s, u, v
integer, intent(in) :: pot
real(8), intent(out) :: ex, ec, vx, vc
! constants
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: g = 0.8098d0
real(8), parameter :: a = 2.804d0
real(8), parameter :: c = 0.7168d0
! local variables
real(8) s2,exlda, vxlda, eclda, vclda, X, Xs, Xss
real(8) F, Fs, Fss, Hx, Hxs, Hxss, Hc, Hcs, Hcss
real(8) zb, zbs, zbss, w
real(8) n0b, n0bs, n0bss
real(8) ln0b, ln0bs, ln0bss
real(8) zbb, zbbc, zbbs, zbbss
real(8) fxb, fxbs, fxbss
! cutoff
if((rho .le. 1.0d-16)) then
   ex = 0.0d0
   ec = 0.0d0
   vx = 0.0d0
   vc = 0.0d0
   return
endif
s2 = s**2
! LDA correlation
call xc_am05_ldapwc(rho,eclda,vclda)
! LDA exchange
call xc_am05_ldax(rho,exlda,vxlda)
!------------------!
!     exchange     !
!------------------!
! interpolation index
X = 1.0d0 - a*s2/(1.0d0 + a*s2)
! Airy LAA refinement function
call xc_am05_labertw(s**(3.0d0/2.0d0)/sqrt(24.0d0),w)
zb = (3.0d0/2.0d0*w)**(2.0d0/3.0d0)
n0b = w/(2.0d0*pi**2*s**3)
ln0b = -3.0d0/(2.0d0*pi)*(3.0d0*pi**2*n0b)**(1.0d0/3.0d0)
zbbc = ((4.0d0/3.0d0)**(1.0d0/3.0d0)*2.0d0*pi/3.0d0)**4
zbb = (zbbc*zb**2 + zb**4)**(1.0d0/4.0d0)
Fxb = -1.0d0/(ln0b*2.0d0*zbb)
F = (c*s2 + 1.0d0)/(c*s2/Fxb + 1.0d0)
! exchange refinement function
Hx = X + (1.0d0 - X)*F
! exchange energy per particle, Ex = Integrate[n*ex]
ex = exlda*Hx
!---------------------!
!     correlation     !
!---------------------!
! correlation refinement function
Hc = X + g*(1.0d0 - X)
! correlation energy per particle, Ec = Integrate[rho*ec]
ec = eclda*Hc
if (pot .eq. 0) return
!----------------------------!
!     exchange potential     !
!----------------------------!
! interpolation index derivatives, dX/ds
Xs = -2.0d0*a*s/(1.0d0 + a*s2)**2
Xss = 2.0d0*a*(3.0d0*a*s2-1.0d0)/(1.0d0+a*s2)**3
! airy LAA refinement function derivatives, dF/ds
zbs = zb/(s + s*w)
zbss = - zb*w*(5.0d0+2.0d0*w)/(2.0d0*s2*(1.0d0+w)**3)
n0bs = sqrt(zb)*(-2.0d0*zb+s*zbs)/(2.0d0*pi**2*s2**2)
n0bss = (16.0d0*zb**2+s**2*zbs**2+2.0d0*s*zb*(-6.0d0* &
     zbs+s*zbss))/(4.0d0*pi**2*s**5*sqrt(zb))
ln0bs = -(3.0d0/pi)**(1.0d0/3.0d0)*n0bs/ &
     (2.0d0*n0b**(2.0d0/3.0d0))
ln0bss = (2.0d0*n0bs**2-3.0d0*n0b*n0bss)/(2.0d0* &
     3.0d0**(2.0d0/3.0d0)*pi**(1.0d0/3.0d0)*n0b**(5.0d0/3.0d0))
zbbs = zb*(zbbc+2*zb**2)*zbs/ &
     (2.0d0*(zb**2*(zbbc+zb**2))**(3.0d0/4.0d0))
zbbss = zb**2*(-zbbc*(zbbc-2.0d0*zb**2)*zbs**2+ &
     2.0d0*zb*(zbbc+zb**2)*(zbbc+2.0d0*zb**2)*zbss)/ &
     (4.0d0*(zb**2*(zbbc+zb**2))**(7.0d0/4.0d0))
Fxbs = (zbb*ln0bs+ln0b*zbbs)/(2.0d0*ln0b**2*zbb**2)
Fxbss = (-2.0d0*ln0b**2*zbbs**2+zbb**2*(-2.0d0*ln0bs**2 + &
     ln0b*ln0bss)+ln0b*zbb*(-2.0d0*ln0bs*zbbs+ln0b*zbbss))/ &
     (2.0d0*ln0b**3*zbb**3)
Fs = (c*s*(2.0d0*(Fxb-1.0d0)*Fxb + s*(1.0d0+c*s2)*Fxbs))/ &
     (c*s2 + Fxb)**2
Fss = (c*(-2.0d0*(3.0d0*c*s2-Fxb)*(Fxb-1.0d0)*Fxb+ &
     4.0d0*s*(-c*s2+Fxb+2.0d0*c*s2*Fxb)*Fxbs - &
     2.0d0*s2*(1.0d0+c*s2)*Fxbs**2+s2*(1.0d0+c*s2)* &
     (c*s2 + Fxb)*Fxbss))/(c*s2+Fxb)**3
! GGA refinement function derivatives, dF/ds
Hxs = - (X - 1.0d0)*Fs - (F - 1.0d0)*Xs
Hxss = - 2.0d0*Fs*Xs - (X - 1.0d0)*Fss - (F - 1.0d0)*Xss
! vx formula for gradient dependent functional,
! generalized form of Eq. (24) in PRB 33, 8800 (1986)
vx = vxlda*(Hx - s*Hxs) + &
     exlda*((4.0d0/3.0d0*s-v/s)*Hxs - &
     (u-4.0d0/3.0d0*s**3)*(Hxss/s-Hxs/s2))
!-------------------------------!
!     correlation potential     !
!-------------------------------!
! correlation refinement function derivatives, dF/ds
Hcs = Xs - g*Xs
Hcss = Xss - g*Xss
! vc formula for gradient dependent functional,
! generalized form of Eq. (24) in Phys. Rev. B 33, 8800 (1986)
vc = vclda*(Hc - s*Hcs) + &
     eclda*((4.0d0/3.0d0*s - v/s)*Hcs - &
     (u - 4.0d0/3.0d0*s**3)*(Hcss/s - Hcs/s2))
return
end subroutine
!EOC

!BOP
! !ROUTINE: xc_am05_ldax
! !INTERFACE:
subroutine xc_am05_ldax(n,ex,vx)
! !INPUT/OUTPUT PARAMETERS:
!   n  : electron density (in,real)
!   ex : exchange energy per electron (out,real)
!   vx : exchange potential (out,real)
! !DESCRIPTION:
!   Local density approximation exchange.
!
! !REVISION HISTORY:
!   Created April 2005 (RAR)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: n
real(8), intent(out) :: ex
real(8), intent(out) :: vx
! constants
real(8), parameter :: pi=3.1415926535897932385d0
vx=-(3.d0*n/pi)**(1.d0/3.d0)
ex=(3.d0/4.d0)*vx
return
end subroutine
!EOC

!BOP
! !ROUTINE: xc_am05_ldapwc
! !INTERFACE:
subroutine xc_am05_ldapwc(n,ec,vc)
! !INPUT/OUTPUT PARAMETERS:
!   n  : electron density (in,real)
!   ec : correlation energy per electron (out,real)
!   vc : correlation potential (out,real)
! !DESCRIPTION:
!   Correlation energy and potential of the Perdew-Wang parameterisation of
!   the Ceperley-Alder electron gas {\it Phys. Rev. B} {\bf 45}, 13244 (1992)
!   and {\it Phys. Rev. Lett.} {\bf 45}, 566 (1980). This is a clean-room
!   implementation from paper.
!
! !REVISION HISTORY:
!   Created April 2005 (RAR)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: n
real(8), intent(out) :: ec
real(8), intent(out) :: vc
! constants
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: a01 = 0.21370d0
real(8), parameter :: b01 = 7.5957d0
real(8), parameter :: b02 = 3.5876d0
real(8), parameter :: b03 = 1.6382d0
real(8), parameter :: b04 = 0.49294d0
! paper actually use this:
! real(8), parameter (A0 = 0.031091d0)
! but routines now "defacto standard" was distributed using:
real(8), parameter :: A0 = 0.0310907d0
! local variables
real(8) rsq
real(8) Q0, Q1, Q1p, ecrs
rsq = (3.0d0/(4.0d0*pi*n))**(1.0d0/6.0d0)
ec = -2.0d0*A0*(1.0d0 + a01*rsq**2)* &
     log(1.0d0 + 1.0d0/ &
          (2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))))
Q0 = -2.0d0*A0*(1.0d0 + a01*rsq**2)
Q1 = 2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))
Q1p = A0*(b01/rsq+2.0d0*b02+3.0d0*b03*rsq+4.0d0*b04*rsq**2)
ecrs = -2.0d0*A0*a01*log(1.0d0 + 1.0d0/Q1)-Q0*Q1p/(Q1**2+Q1)
vc = ec - rsq**2/3.0d0*ecrs
end subroutine
!EOC

!BOP
! !ROUTINE: xc_am05_labertw
! !INTERFACE:
subroutine xc_am05_labertw(z,val)
! !INPUT/OUTPUT PARAMETERS:
!   z   : function argument (in,real)
!   val : value of lambert W function of z (out,real)
! !DESCRIPTION:
!   Lambert $W$-function using the method of Corless, Gonnet, Hare, Jeffrey and
!   Knuth, {\it Adv. Comp. Math.} {\bf 5}, 329 (1996). The approach is based
!   loosely on that in GNU Octave by N. N. Schraudolph, but this implementation
!   is for real values and the principal branch only.
!
! !REVISION HISTORY:
!   Created April 2005 (RAR)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: z
real(8), intent(out) :: val
! local variables
real(8) e,t,p
integer i
! if z too low, go with the first term of the power expansion, z
if (z.lt.1.d-20) then
  val=z
  return
end if
e=exp(1.d0)
! inital guess
if (abs(z+1.d0/e).gt.1.45d0) then
! asymptotic expansion at 0 and Inf
  val=log(z)
  val=val-log(val)
else
! series expansion about -1/e to first order
  val=1.d0*sqrt(2.d0*e*z+2.d0)-1.d0
end if
! find val through iteration
do i=1,10
  p=exp(val)
  t=val*p-z
  if (val.ne.-1.d0) then
    t=t/(p*(val+1.d0)-0.5d0*(val+2.d0)*t/(val+1.d0))
   else
    t=0.d0
   end if
   val=val-t
   if (abs(t).lt.(2.48d0*1.d-14)*(1.d0+abs(val))) return
end do
! this should never happen!
write(*,*)
write(*,'("Error(xc_am05_labertw): iteration limit reached")')
write(*,'(" Likely cause: improper numbers (INFs, NaNs) in density")')
write(*,*)
stop
end subroutine
!EOC

