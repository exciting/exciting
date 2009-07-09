



subroutine hamiltonoverlapocopy_UL(system)
use  modfvsystem
implicit none
type (evsystem) system
integer ::i, ohrank
complex(8), allocatable::tmp(:)
complex(8), pointer::h(:, :), o(:, :)
if(ispacked(system%hamilton).eqv..false.)then
allocate(tmp(getrank(system%hamilton)))
ohrank=getrank(system%hamilton)
h=>system%hamilton%za
o=>system%overlap%za
do i=1, ohrank-1
call zcopy(ohrank-i, h(i, i+1), ohrank, tmp, 1)
tmp(1:ohrank-1)=conjg(tmp(1:ohrank-1))
call zcopy(ohrank-i, tmp, 1, h(i+1, i), 1)
end do

do i=1, ohrank-1
call zcopy(ohrank-i, o(i, i+1), ohrank, tmp, 1)
tmp(1:ohrank-1)=conjg(tmp(1:ohrank-1))
call zcopy(ohrank-i, tmp, 1, o(i+1, i), 1)
end do
endif
end subroutine 
