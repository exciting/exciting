
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rschrod
! !INTERFACE:
subroutine rschrod(n,l,np,nr,r,vr,eval,p0,q0)
! !INPUT/OUTPUT PARAMETERS:
!   n    : principal quantum number (in,integer)
!   l    : quantum number l (in,integer)
!   np   : order of predictor-corrector polynomial (in,integer)
!   nr   : number of radial mesh points (in,integer)
!   r    : radial mesh (in,real(nr))
!   vr   : potential on radial mesh (in,real(nr))
!   eval : eigenvalue without rest-mass energy (inout,real)
!   p0   : P component of the radial wavefunction (out,real(nr))
!   q0   : Q component of the radial wavefunction (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the scalar relativistic radial Schr\"{o}dinger
!   equation for a given potential $v(r)$ and quantum numbers $n$ and $l$. The
!   method involves integrating the equation using the predictor-corrector
!   method and adjusting $E$ until the number of nodes in the wavefunction
!   equals $n-l-1$. The calling routine must provide an initial estimate for the
!   eigenvalue. Note that the arrays {\tt P0} and {\tt Q0} represent the radial
!   functions multiplied by $r$. See routine {\tt rschrodint}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
integer, intent(in) :: l
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
real(8), intent(inout) :: eval
real(8), intent(out) :: p0(nr)
real(8), intent(out) :: q0(nr)
! local variables
integer, parameter :: maxit=2000
integer it,nn,ir,irm,nnd,nndp
! energy convergence tolerance
real(8), parameter :: eps=1.d-14
real(8) t1,de
! automatic arrays
real(8) p1(nr),q1(nr),fr(nr),gr(nr),cf(3,nr)
if (n.lt.1) then
  write(*,*)
  write(*,'("Error(rschrod): n < 1 : ",I8)') n
  write(*,*)
  stop
end if
if (l.gt.n-1) then
  write(*,*)
  write(*,'("Error(rschrod): l > n-1 : ",2I8)') l,n
  write(*,*)
  stop
end if
de=1.d0
nndp=0
do it=1,maxit
! integrate the Schrodinger equation
  call rschroddme(0,l,eval,np,nr,r,vr,nn,p0,p1,q0,q1)
! check the number of nodes
  nnd=nn-(n-l-1)
  if (nnd.gt.0) then
    eval=eval-de
  else
    eval=eval+de
  end if
  if (it.gt.1) then
    if ((nnd.ne.0).or.(nndp.ne.0)) then
      if (nnd*nndp.le.0) de=de*0.5d0
    end if
  end if
  nndp=nnd
  if (de.lt.eps*(abs(eval)+1.d0)) goto 10
end do
write(*,*)
write(*,'("Error(rschrod): maximum iterations exceeded")')
write(*,*)
stop
10 continue
! find effective infinity and set wavefunction to zero after that point
irm=nr
do ir=2,nr
  if ((p0(ir-1)*p0(ir).lt.0.d0).or.(p1(ir-1)*p1(ir).lt.0.d0)) irm=ir
end do
p0(irm:nr)=0.d0
! normalise (major component only)
do ir=1,nr
  fr(ir)=p0(ir)**2
end do
call fderiv(-1,nr,r,fr,gr,cf)
t1=sqrt(abs(gr(nr)))
if (t1.gt.0.d0) then
  t1=1.d0/t1
else
  write(*,*)
  write(*,'("Error(rschrod): zero wavefunction")')
  write(*,*)
  stop
end if
p0(:)=t1*p0(:)
q0(:)=t1*q0(:)
return
end subroutine
