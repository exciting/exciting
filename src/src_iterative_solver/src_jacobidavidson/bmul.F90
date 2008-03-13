subroutine bmul(n,q,r)
use jacobidavidsoncommon
 use modfvsystem
use modmain,only:zzero,zone
implicit none
integer ,intent(in)::n
complex(8),intent(out)::r(n)
complex(8),intent(in)::q(n)

	call	Hermiteanmatrixvector(system%overlap,zone,q,zzero,r)
	
end subroutine bmul
