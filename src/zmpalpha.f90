
subroutine hupdate(z,i,j)
use modfvsystem
integer ,intent(in)::i,j
complex(8),intent(in)::z
integer::ipx
if(packed)then
ipx=((i-1)*i)/2 + j
hp(ipx)=hp(ipx)+z
else
h(i,j)=h(i,j)+z
endif
return
end subroutine

subroutine oupdate(z,i,j)
use modfvsystem
integer ,intent(in)::i,j
complex(8),intent(in)::z
if(packed)then
ipx=((i-1)*i)/2 + j
op(ipx)=op(ipx)+z
else
o(i,j)=h(i,j)+z
endif

return
end subroutine
