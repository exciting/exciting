subroutine bmul(n,q,r)
use jacobidavidsoncommon
 use modfvsystem
use modmain,only:zzero,zone
implicit none
integer ,intent(in)::n
complex(8),intent(in)::r(n)
complex(8),intent(out)::q(n)
	call zhpmv("U",n,zone,overlapp,q(1), 1,zzero,r(1), 1)
	
end subroutine bmul
