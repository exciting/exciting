
subroutine zmpalpha(a,z,i,j)
integer ,intent(in)::i,j
complex(8),intent(in)::z
complex(8),intent(inout)::a(*)
integer::ipx
ipx=((i-1)*i)/2 + j
a(ipx)=a(ipx)+z
return
end subroutine

subroutine zmpalphanp(a,z,i,j)
integer ,intent(in)::i,j
complex(8),intent(in)::z
complex(8),intent(inout)::a(:,:)
a(i,j)=a(i,j)+z
return
end subroutine
