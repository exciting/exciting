subroutine rayleighqotient(n,v,h,o,e)
integer, intent(in)::n
complex(8) ,intent(in)::h(n*(n+1)/2), o(n*(n+1)/2),v(n)
real ,intent(out)::e
complex(8) zdotc
external zdotc
complex(8) :: vwork(n)
complex(8)::vhv,vov
call zhpmv("U",n,dcmplx(1.0,0.0),h,v,1, dcmplx(0,0),vwork,1)
vhv= zdotc(n,v,1,vwork,1)
call zhpmv("U",n,dcmplx(1.0,0.0),o,v,1, dcmplx(0,0),vwork,1)
vhv= zdotc(n,v,1,vwork,1)
e=vhv/vov
end subroutine
