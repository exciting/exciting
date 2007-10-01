
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rschrodint
! !INTERFACE:
subroutine rschrodint(m,l,e,np,nr,r,vr,nn,p0p,p0,p1,q0,q1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   l   : angular momentum quantum number (in,integer)
!   e   : energy (in,real)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   p0p : m-1 th energy derivative of P (in,real(nr))
!   p0  : m th energy derivative of P (out,real(nr))
!   p1  : radial derivative of p0 (out,real(nr))
!   q0  : m th energy derivative of Q (out,real(nr))
!   q1  : radial derivative of q0 (out,real(nr))
! !DESCRIPTION:
!   Integrates the $m$th energy derivative of the scalar relativistic radial
!   Schr\"{o}dinger equation from $r=0$ outwards. This involves using the
!   predictor-corrector method to solve the coupled first-order equations (in
!   atomic units)
!   \begin{align*}
!    \frac{d}{dr}P^{(m)}_l&=2MQ^{(m)}_l+\frac{1}{r}P^{(m)}_l\\
!    \frac{d}{dr}Q^{(m)}_l&=-\frac{1}{r}Q^{(m)}_l+\left[\frac{l(l+1)}{2Mr^2}
!    +(V-E)\right]P^{(m)}_l-mP^{(m-1)}_l,
!   \end{align*}
!   where $V$ is the external potential, $E$ is the eigenenergy and
!   $M=1-V/2c^2$. Following the convention of Koelling and Harmon, {\it J. Phys.
!   C: Solid State Phys.} {\bf 10}, 3107 (1977), the functions $P_l$ and $Q_l$
!   are defined by
!   \begin{align*}
!    P_l&=rg_l\\
!    Q_l&=\frac{r}{2M}\frac{dg_l}{dr},
!   \end{align*}
!   where $g_l$ is the major component of the Dirac equation (see the routine
!   {\tt rdiracint}). Note that in order to make the equations linear in energy,
!   the full definition $M=1+(E-V)/2c^2$ is not used. If $m=0$ then the array
!   {\tt p0p} is not referenced.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m
integer, intent(in) :: l
real(8), intent(in) :: e
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
integer, intent(out) :: nn
real(8), intent(in) :: p0p(nr)
real(8), intent(out) :: p0(nr)
real(8), intent(out) :: p1(nr)
real(8), intent(out) :: q0(nr)
real(8), intent(out) :: q1(nr)
! local variables
integer ir,ir0,npl
! fine-structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
real(8) rm,ri,t1
! automatic arrays
real(8) c(np)
! external functions
real(8) polynom
external polynom
! estimate r -> 0 boundary values
q0(1)=0.d0
p1(1)=1.d0
q1(1)=0.d0
rm=1.d0-0.5d0*(alpha**2)*vr(1)
t1=dble(l*(l+1))/(2.d0*rm*r(1)**2)
p0(1)=r(1)*(p1(1)-2.d0*rm*q0(1))
q0(1)=r(1)*((t1+vr(1)-e)*p0(1)-q1(1))
if (m.ne.0) then
  q1(1)=q1(1)-dble(m)*p0p(1)
end if
nn=0
do ir=2,nr
  rm=1.d0-0.5d0*(alpha**2)*vr(ir)
  ri=1.d0/r(ir)
  t1=dble(l*(l+1))/(2.d0*rm*r(ir)**2)
! predictor-corrector order
  npl=min(ir,np)
  ir0=ir-npl+1
  p1(ir)=polynom(0,npl-1,r(ir0),p1(ir0),c,r(ir))
  q1(ir)=polynom(0,npl-1,r(ir0),q1(ir0),c,r(ir))
! integrate to find wavefunction
  p0(ir)=polynom(-1,npl,r(ir0),p1(ir0),c,r(ir))+p0(ir0)
  q0(ir)=polynom(-1,npl,r(ir0),q1(ir0),c,r(ir))+q0(ir0)
! compute the derivatives
  p1(ir)=2.d0*rm*q0(ir)+p0(ir)*ri
  q1(ir)=(t1+vr(ir)-e)*p0(ir)-q0(ir)*ri
  if (m.ne.0) then
    q1(ir)=q1(ir)-dble(m)*p0p(ir)
  end if
! integrate for correction
  p0(ir)=polynom(-1,npl,r(ir0),p1(ir0),c,r(ir))+p0(ir0)
  q0(ir)=polynom(-1,npl,r(ir0),q1(ir0),c,r(ir))+q0(ir0)
! compute the derivatives again
  p1(ir)=2.d0*rm*q0(ir)+p0(ir)*ri
  q1(ir)=(t1+vr(ir)-e)*p0(ir)-q0(ir)*ri
  if (m.ne.0) then
    q1(ir)=q1(ir)-dble(m)*p0p(ir)
  end if
! check for overflow
  if ((abs(p0(ir)).gt.1.d100).or.(abs(p1(ir)).gt.1.d100).or. &
      (abs(q0(ir)).gt.1.d100).or.(abs(q1(ir)).gt.1.d100)) then
    p0(ir:nr)=p0(ir)
    p1(ir:nr)=p1(ir)
    q0(ir:nr)=q0(ir)
    q1(ir:nr)=q1(ir)
    return
  end if
! check for node
  if (p0(ir-1)*p0(ir).lt.0.d0) nn=nn+1
end do
return
end subroutine
!EOC
