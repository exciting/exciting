subroutine amul(n,q,r)
use jacobidavidsoncommon
use modmain,only:zzero,zone
implicit none
integer ,intent(in)::n
complex(8),intent(out)::r(n)
complex(8),intent(in)::q(n)
call	Hermiteanmatrixvector(system%hamilton,zone,q,zzero,r)
end subroutine amul
