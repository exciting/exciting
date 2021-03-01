
! Copyright (C) 2007 J. K. Dewhurst and D. W. H. Rankin.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine minf_nm(n,x,maxit,iter,eps)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(inout) :: x(n,n+1)
integer, intent(in) :: maxit
integer, intent(out) :: iter
real(8), intent(in) :: eps
! local variables
integer i,j,il,iu
! Nelder-Mead parmeters
real(8), parameter :: alpha=1.d0
real(8), parameter :: gamma=2.d0
real(8), parameter :: beta=0.5d0
real(8), parameter :: sigma=0.5d0
real(8) fr,fe,fc,sum,t1
! automatic arrays
real(8) f(n+1),xm(n),xr(n),xe(n),xc(n)
! external functions
real(8) fmin_nm
external fmin_nm
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(minf_nm): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
! evaluate the function at each vertex
do i=1,n+1
  f(i)=fmin_nm(x(1,i))
end do
iter=0
10 continue
iter=iter+1
if (iter.ge.maxit) return
! find the lowest and highest vertex
il=1
iu=1
do i=2,n+1
  if (f(i).lt.f(il)) il=i
  if (f(i).gt.f(iu)) iu=i
end do
! check for convergence
if ((f(iu)-f(il)).lt.eps) return
! compute the mean of the n lowest vertices
t1=1.d0/dble(n)
do i=1,n
  sum=0.d0
  do j=1,iu-1
    sum=sum+x(i,j)
  end do
  do j=iu+1,n+1
    sum=sum+x(i,j)
  end do
  xm(i)=t1*sum
end do
xr(:)=xm(:)+alpha*(xm(:)-x(:,iu))
fr=fmin_nm(xr)
if (f(il).gt.fr) goto 30
if ((f(il).le.fr).and.(fr.lt.f(iu))) then
! reflection
  x(:,iu)=xr(:)
  f(iu)=fr
  goto 10
else
  goto 40
end if
30 continue
xe(:)=xm(:)+gamma*(xr(:)-xm(:))
fe=fmin_nm(xe)
if (fr.gt.fe) then
! expansion
  x(:,iu)=xe(:)
  f(iu)=fe
else
! reflection
  x(:,iu)=xr(:)
  f(iu)=fr
end if
goto 10
40 continue
xc(:)=xm(:)+beta*(x(:,iu)-xm(:))
fc=fmin_nm(xc)
if (fc.lt.f(iu)) then
! contraction
  x(:,iu)=xc(:)
  f(iu)=fc
  goto 10
end if
! shrinkage
do j=1,il-1
  x(:,j)=x(:,il)+sigma*(x(:,j)-x(:,il))
  f(j)=fmin_nm(x(1,j))
end do
do j=il+1,n+1
  x(:,j)=x(:,il)+sigma*(x(:,j)-x(:,il))
  f(j)=fmin_nm(x(1,j))
end do
goto 10
return
end subroutine

