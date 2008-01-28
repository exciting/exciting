subroutine initprecon(n,sigma,h,o)
use modmain,only:zzero,zone
use jacobidavidsoncommon
implicit none
integer,intent(in)::n
complex(8),intent(in)::sigma
complex(8),intent(in):: h(n*(n+1)/2),o(n*(n+1)/2)

integer np,info
np=n*(n+1)/2
 allocate(p(np),ipiv(n))
  call zcopy(np,h,1,p,1)
  call zaxpy(np,-sigma,o,1,p,1)
  call zhptrf('U', n, p, IPIV, info )
end subroutine

subroutine deallocprecon()
use jacobidavidsoncommon
implicit none
deallocate(p,ipiv)
end subroutine