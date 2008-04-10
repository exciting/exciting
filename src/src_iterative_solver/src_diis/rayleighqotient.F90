subroutine rayleighqotient(n,m,evecfv, h,s,evalfv)
 
  implicit none
  integer, intent(in)::n,m
  complex(8) ,intent(in)::h(n,m),s(n,m),evecfv(n,m)
  real(8) ,intent(out)::evalfv(m)
  complex(8) zdotc
  external zdotc
  complex(8) :: vwork(n)
  complex(8)::vhv,vsv
  integer i
  do i=1,m
     vhv= zdotc(n,evecfv(1,i),1,h(1,i),1)
     vsv= zdotc(n,evecfv(1,i),1,s(1,i),1)
     evalfv(i)=vhv/vsv
	 write (*,*) " vhv ,shv",vhv ,vsv
  end do
end subroutine rayleighqotient
