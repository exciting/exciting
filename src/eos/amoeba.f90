subroutine amoeba(ftol,nparam,p,iter)
implicit none
! arguments
real(8), intent(in) :: ftol
integer, intent(in) :: nparam
real(8), intent(inout) :: p(nparam,nparam+1)
integer, intent(out) :: iter
! local variables
integer, parameter :: maxit=1000000
integer i,ihi,ilo,inhi,j,m,n
real(8) rtol,rsum,swap,ysave,ytry
! external functions
real(8) amotry,fopt
external amotry,fopt
! automatic arrays
real(8) psum(nparam),y(nparam+1)
do i=1,nparam+1
  y(i)=fopt(p(1,i))
end do
iter=0
10 continue
do n=1,nparam
  rsum=0.d0
  do m=1,nparam+1
    rsum=rsum+p(n,m)
  end do
  psum(n)=rsum
end do
20 continue
ilo=1
if (y(1).gt.y(2)) then
  ihi=1
  inhi=2
else
  ihi=2
  inhi=1
end if
do i=1,nparam+1
  if (y(i).le.y(ilo)) ilo=i
  if (y(i).gt.y(ihi)) then
  inhi=ihi
  ihi=i
  else if (y(i).gt.y(inhi)) then
    if (i.ne.ihi) inhi=i
  end if
end do
rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
if (rtol.lt.ftol) then
  swap=y(1)
  y(1)=y(ilo)
  y(ilo)=swap
  do n=1,nparam
    swap=p(n,1)
    p(n,1)=p(n,ilo)
    p(n,ilo)=swap
  end do
  return
end if
if (iter.ge.maxit) then
  write(*,*)
  write(*,'("Error(amoeba): maxit exceeded : ",I10)') maxit
  write(*,*)
  stop
end if
iter=iter+2
ytry=amotry(p,y,psum,nparam,ihi,-1.d0)
if (ytry.le.y(ilo)) then
  ytry=amotry(p,y,psum,nparam,ihi,2.d0)
else if (ytry.ge.y(inhi)) then
  ysave=y(ihi)
  ytry=amotry(p,y,psum,nparam,ihi,0.5d0)
  if (ytry.ge.ysave) then
    do i=1,nparam+1
      if (i.ne.ilo) then
        do j=1,nparam
          psum(j)=0.5d0*(p(j,i)+p(j,ilo))
          p(j,i)=psum(j)
        end do
        y(i)=fopt(psum)
      end if
    end do
    iter=iter+nparam
    goto 10
  end if
else
  iter=iter-1
end if
goto 20
end subroutine
