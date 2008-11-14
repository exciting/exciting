
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixadapt
! !INTERFACE:
subroutine mixadapt(iscl,beta0,betainc,betadec,n,nu,mu,beta,f,d)
! !INPUT/OUTPUT PARAMETERS:
!   iscl    : self-consistent loop number (in,integer)
!   beta0   : initial value for mixing parameter (in,real)
!   betainc : mixing parameter increase (in,real)
!   betadec : mixing parameter decrease (in,real)
!   n       : vector length (in,integer)
!   nu      : current output vector as well as the next input vector of the
!             self-consistent loop (inout,real(n))
!   mu      : used internally (inout,real(n))
!   beta    : used internally (inout,real(n))
!   f       : used internally (inout,real(n))
!   d       : RMS difference between old and new output vectors (out,real)
! !DESCRIPTION:
!   Given the input vector $\mu^i$ and output vector $\nu^i$ of the $i$th
!   self-consistent loop, this routine generates the next input vector to the
!   loop using an adaptive mixing scheme. The $j$th component of the output
!   vector is mixed with a fraction of the same component of the input vector:
!   $$ \mu^{i+1}_j=\beta^i_j\nu^i_j+(1-\beta^i_j)\mu^i_j, $$
!   where $\beta^i_j$ is set to $\beta_0$ at initialisation and increased by
!   scaling with $\beta_{\rm inc}$ ($>1$) if $f^i_j\equiv\nu^i_j-\mu^i_j$ does
!   not change sign between loops. If $f^i_j$ does change sign, then $\beta^i_j$
!   is scaled by $\beta_{\rm dec}$ ($<1$). Note that the array {\tt nu} serves
!   for both input and output, and the arrays {\tt mu}, {\tt beta} and {\tt f}
!   are used internally and should not be changed between calls. The routine is
!   initialised at the first iteration and is thread-safe so long as each thread
!   has its own independent work array. Complex arrays may be passed as real
!   arrays with $n$ doubled.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Modified, September 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: iscl
real(8), intent(in) :: beta0
real(8), intent(in) :: betainc
real(8), intent(in) :: betadec
integer, intent(in) :: n
real(8), intent(inout) :: nu(n)
real(8), intent(inout) :: mu(n)
real(8), intent(inout) :: beta(n)
real(8), intent(inout) :: f(n)
real(8), intent(out) :: d
! local variables
integer i
real(8) t1
if (iscl.le.1) then
  mu(:)=nu(:)
  f(:)=0.d0
  beta(:)=beta0
  d=1.d0
  return
end if
do i=1,n
  t1=nu(i)-mu(i)
  if (t1*f(i).gt.0.d0) then
    beta(i)=beta(i)*betainc
    if (beta(i).gt.1.d0) beta(i)=1.d0
  else
    beta(i)=beta(i)*betadec
  end if
  f(i)=t1
end do
nu(:)=beta(:)*nu(:)+(1.d0-beta(:))*mu(:)
d=0.d0
do i=1,n
  d=d+f(i)**2
end do
d=sqrt(d/dble(n))
mu(:)=nu(:)
return
end subroutine
!EOC

