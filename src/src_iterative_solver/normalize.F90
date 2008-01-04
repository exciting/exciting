subroutine  normalize(n,m,overlap,evecfv)	
use modmain,only:zone,zzero
use diisinterfaces
implicit none
integer, intent (in)::n,m
complex(8),intent(in)::overlap(n,n)
complex(8),intent(out)::evecfv(n,m)
complex(8)::tmp(n),z
integer i
do i=1,m
		call zhemv('U',n,zone,overlap(1,1),n,evecfv(1,i),1,&
      	zzero,tmp,1)
      	z =sqrt(dble (zdotc(n,evecfv(1,i),1,tmp,1)))
      	z=1.0/z
      	call zscal(n,z,evecfv(1,i),1)
end do
end subroutine
