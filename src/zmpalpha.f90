
subroutine hupdate(z,i,j)
use modfvsystem,only:hamilton,hamiltonp,packed
integer ,intent(in)::i,j
complex(8),intent(in)::z
integer::ipx
if(packed)then
ipx=((i-1)*i)/2 + j
hamiltonp(ipx)=hamiltonp(ipx)+z
else
hamilton(i,j)=hamilton(i,j)+z
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
overlap(i,j)=overlap(i,j)+z
endif

return
end subroutine
