subroutine hamiltonoverlapocopy_UL
use  modfvsystem,only:hamilton,overlap,ohrank
integer ::i
complex(8)::tmp(ohrank)
do i=1, ohrank-1
call zcopy(ohrank-i,hamilton(i,i+1),ohrank,tmp,1)
call zcopy(ohrank-i,conjg(tmp),1,hamilton(i+1,i),1)
end do

do i=1, ohrank-1
call zcopy(ohrank-i,overlap(i,i+1),ohrank,tmp,1)
call zcopy(ohrank-i,conjg(tmp),1,overlap(i+1,i),1)
end do
end subroutine 