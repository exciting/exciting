subroutine fitdata
use modeos
implicit none
! local variables
integer iter,ip
real(8) ftol
! automatic arrays
real(8) p(nparam,nparam+1)
! external functions
real(8) fopt
external fopt
ftol=1.d-9
iter=0
! initial guess: it is assumed that param(1)=V0, param(2)=E0 and param(3)=B0
p(:,1)=0.d0
p(1,1)=vpt(1)
p(2,1)=ept(1)
p(3,1)=0.003d0
! fit V0 and E0
do ip=1,nparam
  p(:,ip+1)=p(:,1)
end do
p(1,2)=p(1,2)+1.d0
p(2,3)=p(2,3)+0.1d0
call amoeba(ftol,nparam,p,iter)
! fit V0, E0 and B0
do ip=1,nparam
  p(:,ip+1)=p(:,1)
end do
p(1,2)=p(1,2)+1.d0
p(2,3)=p(2,3)+0.1d0
p(3,4)=p(3,4)+0.001d0
call amoeba(ftol,nparam,p,iter)
! fit everything
do ip=1,nparam
  p(:,ip+1)=p(:,1)
  p(ip,ip+1)=p(ip,ip+1)+0.1d0
end do
call amoeba(ftol,nparam,p,iter)
popt(1:nparam)=p(1:nparam,1)
return
end subroutine
