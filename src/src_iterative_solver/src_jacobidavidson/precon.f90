subroutine precon(n,q)
use jacobidavidsoncommon
use modmain,only:zzero,zone
implicit none
integer ,intent(in)::n
complex(8),intent(inout)::q(n)
integer info
call zhptrs('U', n, 1, p(1), IPIV, q(1), n, INFO )
end subroutine
