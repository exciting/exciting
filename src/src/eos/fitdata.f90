
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fitdata
use modmain
implicit none
! local variables
integer, parameter :: maxit=1000000
integer i,iter
real(8), parameter :: eps=1.d-14
! automatic arrays
real(8) x(nparam,nparam+1)
! initial guess: it is assumed that param(1)=V0, param(2)=E0 and param(3)=B0
x(:,1)=0.d0
x(1,1)=vpt(1)
x(2,1)=ept(1)
x(3,1)=0.003d0
! fit V0 and E0
do i=1,nparam
  x(:,i+1)=x(:,1)
end do
x(1,2)=x(1,2)+1.d0
x(2,3)=x(2,3)+0.1d0
call minf_nm(nparam,x,maxit,iter,eps)
! fit V0, E0 and B0
do i=1,nparam
  x(:,i+1)=x(:,1)
end do
x(1,2)=x(1,2)+1.d0
x(2,3)=x(2,3)+0.1d0
x(3,4)=x(3,4)+0.001d0
call minf_nm(nparam,x,maxit,iter,eps)
! fit everything
do i=1,nparam
  x(:,i+1)=x(:,1)
  x(i,i+1)=x(i,i+1)+0.1d0
end do
call minf_nm(nparam,x,maxit,iter,eps)
popt(1:nparam)=x(1:nparam,1)
return
end subroutine

