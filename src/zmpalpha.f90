
subroutine hupdate(z,i,j)
use modfvsystem,only:hamilton,hamiltonp,packed
integer ,intent(in)::i,j
complex(8),intent(in)::z
integer::ipx
if(packed)then
ipx=((i-1)*i)/2 + j
hamiltonp(ipx)=hamiltonp(ipx)+z
else
if(j.le.i)then
hamilton(j,i)=hamilton(j,i)+z
else
write(*,*)"warning lower part of hamilton updated"
endif
endif
return
end subroutine

subroutine oupdate(z,i,j)
use modfvsystem,only:overlapp,overlap,packed
integer ,intent(in)::i,j
complex(8),intent(in)::z
if(packed)then
ipx=((i-1)*i)/2 + j
overlapp(ipx)=overlapp(ipx)+z
else
if(j.le.i)then
overlap(j,i)=overlap(j,i)+z
else
write(*,*)"warning lower part of overlap updated"
endif
endif

return
end subroutine
