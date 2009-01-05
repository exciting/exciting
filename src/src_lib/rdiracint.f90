
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdiracint
! !INTERFACE:
subroutine rdiracint(m,kpa,e,np,nr,r,vr,nn,g0p,f0p,g0,g1,f0,f1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   kpa : quantum number kappa (in,integer)
!   e   : energy (in,real)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   g0p : m-1 th energy derivative of the major component multiplied by r
!         (in,real(nr))
!   f0p : m-1 th energy derivative of the minor component multiplied by r
!         (in,real(nr))
!   g0  : m th energy derivative of the major component multiplied by r
!         (out,real(nr))
!   g1  : radial derivative of g0 (out,real(nr))
!   f0  : m th energy derivative of the minor component multiplied by r
!         (out,real(nr))
!   f1  : radial derivative of f0 (out,real(nr))
! !DESCRIPTION:
!   Integrates the $m$th energy derivative of the radial Dirac equation from
!   $r=0$ outwards. This involves using the predictor-corrector method to solve
!   the coupled first-order equations (in atomic units)
!   \begin{align*}
!    \left(\frac{d}{dr}+\frac{\kappa}{r}\right)G^{(m)}_\kappa&=\frac{1}{c}
!    \{2E_0+E-V\}F^{(m)}_\kappa+\frac{m}{c}F^{(m-1)}_\kappa\\
!    \left(\frac{d}{dr}-\frac{\kappa}{r}\right)F^{(m)}_\kappa&=
!    -\frac{1}{c}\{E-V\}G^{(m)}_\kappa-\frac{m}{c}G^{(m-1)}_\kappa,
!   \end{align*}
!   where $G^{(m)}_\kappa=rg^{(m)}_\kappa$ and $F^{(m)}_\kappa=rf^{(m)}_\kappa$
!   are the $m$th energy derivatives of the major and minor components
!   multiplied by $r$, respectively; $V$ is the external potential; $E_0$ is the
!   electron rest energy; $E$ is the eigen energy (excluding $E_0$); and
!   $\kappa=l$ for $j=l-\frac{1}{2}$ or $\kappa=-(l+1)$ for $j=l+\frac{1}{2}$.
!   If $m=0$ then the arrays {\tt g0p} and {\tt f0p} are not referenced.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m
integer, intent(in) :: kpa
real(8), intent(in) :: e
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
integer, intent(out) :: nn
real(8), intent(in) :: g0p(nr)
real(8), intent(in) :: f0p(nr)
real(8), intent(out) :: g0(nr)
real(8), intent(out) :: g1(nr)
real(8), intent(out) :: f0(nr)
real(8), intent(out) :: f1(nr)
! local variables
integer ir,ir0,npl
! fine-structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
! rest mass of electron
real(8), parameter :: e0=1.d0/(alpha**2)
real(8) rkpa,ri
! automatic arrays
real(8) c(np)
! external functions
real(8) polynom
external polynom
if (nr.le.0) then
  write(*,*)
  write(*,'("Error(rdiracint): invalid nr : ",I8)') nr
  write(*,*)
  stop
end if
if ((m.lt.0).or.(m.gt.6)) then
  write(*,*)
  write(*,'("Error(rdiracint): m out of range : ",I8)') m
  write(*,*)
  stop
end if
rkpa=dble(kpa)
! estimate r -> 0 boundary values
f0(1)=0.d0
g1(1)=1.d0
f1(1)=0.d0
ri=1.d0/r(1)
g0(1)=(alpha*(e+2.d0*e0-vr(1))*f0(1)-g1(1))/(rkpa*ri)
f0(1)=(alpha*(e-vr(1))*g0(1)+f1(1))/(rkpa*ri)
if (m.ne.0) then
  g1(1)=g1(1)+alpha*dble(m)*f0p(1)
  f1(1)=f1(1)-alpha*dble(m)*g0p(1)
end if
nn=0
do ir=2,nr
  ri=1.d0/r(ir)
! predictor-corrector order
  npl=min(ir,np)
  ir0=ir-npl+1
  g1(ir)=polynom(0,npl-1,r(ir0),g1(ir0),c,r(ir))
  f1(ir)=polynom(0,npl-1,r(ir0),f1(ir0),c,r(ir))
! integrate to find wavefunction
  g0(ir)=polynom(-1,npl,r(ir0),g1(ir0),c,r(ir))+g0(ir0)
  f0(ir)=polynom(-1,npl,r(ir0),f1(ir0),c,r(ir))+f0(ir0)
! compute the derivatives
  g1(ir)=alpha*(e+2.d0*e0-vr(ir))*f0(ir)-rkpa*ri*g0(ir)
  f1(ir)=-alpha*(e-vr(ir))*g0(ir)+rkpa*ri*f0(ir)
  if (m.ne.0) then
    g1(ir)=g1(ir)+alpha*dble(m)*f0p(ir)
    f1(ir)=f1(ir)-alpha*dble(m)*g0p(ir)
  end if
! integrate for correction
  g0(ir)=polynom(-1,npl,r(ir0),g1(ir0),c,r(ir))+g0(ir0)
  f0(ir)=polynom(-1,npl,r(ir0),f1(ir0),c,r(ir))+f0(ir0)
! compute the derivatives again
  g1(ir)=alpha*(e+2.d0*e0-vr(ir))*f0(ir)-rkpa*ri*g0(ir)
  f1(ir)=-alpha*(e-vr(ir))*g0(ir)+rkpa*ri*f0(ir)
  if (m.ne.0) then
    g1(ir)=g1(ir)+alpha*dble(m)*f0p(ir)
    f1(ir)=f1(ir)-alpha*dble(m)*g0p(ir)
  end if
! check for overflow
  if ((abs(g0(ir)).gt.1.d100).or.(abs(g1(ir)).gt.1.d100).or. &
      (abs(f0(ir)).gt.1.d100).or.(abs(f1(ir)).gt.1.d100)) then
    g0(ir:nr)=g0(ir)
    g1(ir:nr)=g1(ir)
    f0(ir:nr)=f0(ir)
    f1(ir:nr)=f1(ir)
    return
  end if
! check for node
  if (g0(ir-1)*g0(ir).lt.0.d0) nn=nn+1
end do
return
end subroutine
!EOC
