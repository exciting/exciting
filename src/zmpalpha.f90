
subroutine zmpalpha(a,np,z,i,j)
integer ,intent(in)::i,j
complex(8),intent(in)::z
complex(8),intent(inout)::a(np)
integer::ipx
ipx=((i-1)*i)/2 + j
a(ipx)=a(ipx)+z
return
end subroutine

subroutine zmpalphanp(a,n,z,i,j)
integer ,intent(in)::i,j
complex(8),intent(in)::z
complex(8),intent(inout)::a(n,n)
a(i,j)=a(i,j)+z
return
end subroutine
