subroutine rayleighqotient(n,m,evecfv, h,s,evalfv)
use modmain, only: nstfv
implicit none
integer, intent(in)::n,m
complex(8) ,intent(in)::h(n),s(n),evecfv(n,m)
real ,intent(out)::evalfv(m)
complex(8) zdotc
external zdotc
complex(8) :: vwork(n)
complex(8)::vhv,vsv
integer i
do i=1,m
vhv= zdotc(n,evecfv,1,h,1)
vsv= zdotc(n,evecfv,1,s,1)
evalfv(i)=vhv/vsv
end do
end subroutine
