subroutine precon(n,q)
use jacobidavidsoncommon
use modmain,only:zzero,zone
implicit none
integer ,intent(in)::n
complex(8),intent(inout)::q(n)
integer ::i
!call Hermiteanmatrixlinsolve(p,q)
do i=1,n
q(i)=q(i)/pdiag(i)
end do
end subroutine
