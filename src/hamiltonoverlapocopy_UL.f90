subroutine hamiltonoverlapocopy_UL(system)
use  modfvsystem
implicit none
type (evsystem) system
integer ::i,ohrank
complex(8),allocatable::tmp(:)
if(system%hamilton%packed.eqv..false.)then
allocate(tmp(system%hamilton%rank))
ohrank=system%hamilton%rank
do i=1, ohrank-1
call zcopy(ohrank-i,system%hamilton%za(i,i+1),ohrank,tmp,1)
tmp(1:ohrank-1)=conjg(tmp(1:ohrank-1))
call zcopy(ohrank-i,tmp,1,system%hamilton%za(i+1,i),1)
end do

do i=1, ohrank-1
call zcopy(ohrank-i,system%overlap%za(i,i+1),ohrank,tmp,1)
tmp(1:ohrank-1)=conjg(tmp(1:ohrank-1))
call zcopy(ohrank-i,tmp,1,system%overlap%za(i+1,i),1)
end do
endif
end subroutine 